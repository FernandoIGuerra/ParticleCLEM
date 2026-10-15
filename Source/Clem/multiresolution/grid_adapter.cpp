#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include "averager.H"
#include "grid_adapter.H"
#include "line_packing.H"
#include "setup.H"

namespace clem {

namespace {

using ParticleType = ClemContainer::ParticleType;

// the cells of level lev with no line get the clone of the line of their
// parent on level lev - 1. Returns the number of cells filled on this rank
amrex::Long
CloneIntoEmptyCells(ClemContainer& container, const int lev)
{
  auto& pc                       = container.Particles();
  const int crse                 = lev - 1;
  const amrex::BoxArray& fine_ba = pc.ParticleBoxArray(lev);
  const auto& fine_dm            = pc.ParticleDistributionMap(lev);
  const amrex::IntVect ref_ratio = pc.GetParGDB()->refRatio(crse);
  const auto prob_lo             = pc.Geom(lev).ProbLoArray();
  const auto cell_size           = pc.Geom(lev).CellSizeArray();
  const amrex::Real share =
    1.0 / static_cast<amrex::Real>(AMREX_D_TERM(ref_ratio[0], *ref_ratio[1], *ref_ratio[2]));

  // is there any cell to fill
  amrex::Long n_empty = 0;
  for (amrex::MFIter mfi = pc.MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      n_empty += container.ParticlesInCell(lev, mfi.index(), iv).empty() ? 1 : 0;
    }
  }
  amrex::Long n_empty_all = n_empty;
  amrex::ParallelDescriptor::ReduceLongSum(n_empty_all);
  if (n_empty_all == 0) {
    return 0;
  }

  // the lines of the parents, on the rank of their children
  amrex::BoxArray parent_ba = fine_ba;
  parent_ba.coarsen(ref_ratio);
  amrex::MultiFab parents(parent_ba, fine_dm, n_packed, 0);
  parents.setVal(0.0);
  {
    const amrex::MultiFab packed = PackCoarseLines(container, crse, parent_ba);
    parents.ParallelCopy(packed, 0, 0, n_packed, 0, 0, pc.Geom(crse).periodicity());
  }

  std::vector<ParticleType> line;
  std::vector<ParticleType> fresh;
  for (amrex::MFIter mfi = pc.MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    auto& tile           = pc.DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
    const auto& src      = parents.const_array(mfi);

    fresh.clear();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (!container.ParticlesInCell(lev, mfi.index(), iv).empty()) {
        continue;
      }
      const amrex::IntVect iv_crse = amrex::coarsen(iv, ref_ratio);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(src(iv_crse, PackedComp(0, RealData::mass)) > 0.0,
                                       "CLEM regrid: no coarse line under a new fine cell");
      CloneLine(src, iv_crse, share, line);
      for (int k = 0; k < setup::n_lem; ++k) {
        ParticleType& p = line[k];
        p.id()          = ParticleType::NextID();
        p.cpu()         = amrex::ParallelDescriptor::MyProc();
        PlaceOnLine(p, iv, k, setup::n_lem, prob_lo, cell_size);
        fresh.push_back(p);
      }
    }
    // after the loop: push_back may move the storage of the tile
    for (const auto& p : fresh) {
      tile.push_back(p);
    }
  }
  return n_empty;
}

} // namespace

AdaptDiagnostics
AdaptToGrids(ClemContainer& container, const int lbase)
{
  BL_PROFILE("clem::AdaptToGrids()");
  AMREX_ALWAYS_ASSERT(container.IsDefined());

  auto& pc         = container.Particles();
  const int finest = container.FinestLevel();

  AdaptDiagnostics diag;

  // most regrids leave the grids as they were: nothing to do then
  bool changed = static_cast<int>(pc.GetParticles().size()) > finest + 1 &&
                 pc.NumberOfParticlesAtLevel(finest + 1) > 0;
  for (int lev = lbase + 1; lev <= finest; ++lev) {
    changed = changed || container.GridsChanged(lev);
  }
  if (!changed) {
    return diag;
  }

  container.InvalidateCellMaps();
  pc.resizeData(); // a level may be new

  // a level that no longer exists: all its cells left it
  for (int lev = finest + 1; lev < static_cast<int>(pc.GetParticles().size()); ++lev) {
    diag.n_dropped += pc.NumberOfParticlesAtLevel(lev);
    pc.RemoveParticlesAtLevel(lev);
  }

  // from the coarsest changed level up, so that a parent exists before its children
  container.RebuildCellMap(lbase);
  for (int lev = lbase + 1; lev <= finest; ++lev) {
    // a level on unchanged grids, under unchanged parents, is as it was
    if (!container.GridsChanged(lev) && diag.n_refined == 0 && diag.n_dropped == 0) {
      container.RebuildCellMap(lev);
      continue;
    }
    // refile under the new grids; a line whose cell left the level is dropped.
    // Its elements are marked invalid first, so that the redistribute removes
    // them: AMReX aborts on a valid particle outside the grids of its level
    const amrex::Long n_before   = pc.NumberOfParticlesAtLevel(lev);
    const amrex::BoxArray& grids = pc.ParticleBoxArray(lev);
    for (auto& tile : pc.GetParticles(lev)) {
      for (auto& p : tile.second.GetArrayOfStructs()) {
        if (p.id() > 0 && !grids.contains(pc.Index(p, lev))) {
          p.id() = -p.id();
        }
      }
    }
    container.Redistribute(lev);
    diag.n_dropped += n_before - pc.NumberOfParticlesAtLevel(lev);
    container.RebuildCellMap(lev);

    amrex::Long n_refined = CloneIntoEmptyCells(container, lev);
    amrex::ParallelDescriptor::ReduceLongSum(n_refined);
    diag.n_refined += n_refined;
    container.RebuildCellMap(lev);

    const amrex::Long n_expected =
      static_cast<amrex::Long>(setup::n_lem) * pc.ParticleBoxArray(lev).numPts();
    if (pc.NumberOfParticlesAtLevel(lev) != n_expected) {
      amrex::Print() << "CLEM regrid: level " << lev << " holds "
                     << pc.NumberOfParticlesAtLevel(lev) << " elements, expected " << n_expected
                     << '\n';
      amrex::Abort("CLEM regrid: the lines do not fill the new grids");
    }
  }

  // every covered cell of the changed levels as after a sync. A newly covered
  // cell gets the average of its clones, which is the line it held
  for (int lev = finest - 1; lev >= lbase; --lev) {
    AverageDown(container, lev);
  }
  container.RememberGrids();
  return diag;
}

} // namespace clem
