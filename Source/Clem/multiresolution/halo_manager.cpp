#include <algorithm>
#include <cmath>

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_iMultiFab.H>

#include "halo_manager.H"
#include "line_packing.H"
#include "setup.H"

namespace clem {

namespace {

// what a fine cell position around a fine grid is
constexpr int not_halo = 0; // a cell of the fine level, or outside the domain
constexpr int is_halo  = 1; // inside the domain, on the coarser level

} // namespace

void
HaloManager::Clear()
{
  m_lines.clear();
  m_lev = -1;
}

void
HaloManager::FillFromCoarse(ClemContainer& container, const int lev)
{
  BL_PROFILE("clem::HaloManager::FillFromCoarse()");
  AMREX_ALWAYS_ASSERT(container.IsDefined());
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev >= 1 && lev <= container.FinestLevel(),
                                   "FillFromCoarse: lev must be a fine level");
  Clear();

  const int crse                 = lev - 1;
  const auto& pc                 = container.Particles();
  const amrex::Geometry& geom    = pc.Geom(lev);
  const amrex::BoxArray& fine_ba = pc.ParticleBoxArray(lev);
  const auto& fine_dm            = pc.ParticleDistributionMap(lev);
  const amrex::IntVect ref_ratio = pc.GetParGDB()->refRatio(crse);
  // one child of r^d: its share of the extensive part of the coarse line
  const amrex::Real share =
    1.0 / static_cast<amrex::Real>(AMREX_D_TERM(ref_ratio[0], *ref_ratio[1], *ref_ratio[2]));

  if (!container.IsCellMapValid(crse)) {
    container.RebuildCellMap(crse);
  }

  // the coarsened layer of each fine grid: same grid index and rank as the
  // fine grid, so each fine grid reads its own copy
  amrex::BoxArray halo_ba = fine_ba;
  halo_ba.grow(1).coarsen(ref_ratio);

  // 1. and 2.
  amrex::MultiFab halo(halo_ba, fine_dm, n_packed, 0);
  halo.setVal(0.0);
  {
    const amrex::MultiFab packed = PackCoarseLines(container, crse, halo_ba);
    halo.ParallelCopy(packed, 0, 0, n_packed, 0, 0, pc.Geom(crse).periodicity());
  }

  // which positions of the layer are halo cells
  amrex::iMultiFab kind(fine_ba, fine_dm, 1, 1);
  kind.BuildMask(geom.Domain(), geom.periodicity(), not_halo, is_halo, not_halo, not_halo);

  // 3.
  for (amrex::MFIter mfi(halo); mfi.isValid(); ++mfi) {
    const amrex::Box grown_bx = amrex::grow(fine_ba[mfi.index()], 1);
    const auto& src           = halo.const_array(mfi);
    const auto& mask          = kind.const_array(mfi);
    auto& lines               = m_lines[mfi.index()];
    lines.resize(grown_bx);
    lines.setVal(HaloLine{});

    for (amrex::IntVect iv = grown_bx.smallEnd(); iv <= grown_bx.bigEnd(); grown_bx.next(iv)) {
      if (mask(iv) != is_halo) {
        continue;
      }
      const amrex::IntVect iv_crse = amrex::coarsen(iv, ref_ratio);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(src(iv_crse, PackedComp(0, RealData::mass)) > 0.0,
                                       "FillFromCoarse: no coarse line under a halo cell");
      CloneLine(src, iv_crse, share, lines(iv));
    }
  }
  m_lev = lev;
}

bool
HaloManager::IsHaloCell(const int grid_id, const amrex::IntVect& iv) const
{
  const auto it = m_lines.find(grid_id);
  return it != m_lines.end() && it->second.box().contains(iv) && !it->second(iv).empty();
}

const HaloManager::HaloLine&
HaloManager::Line(const int grid_id, const amrex::IntVect& iv) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(IsDefined(), "HaloManager: used before FillFromCoarse");
  return m_lines.at(grid_id)(iv);
}

amrex::Long
HaloManager::NumHaloCells() const
{
  amrex::Long n_cells = 0;
  for (const auto& grid : m_lines) {
    const amrex::Box& bx = grid.second.box();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      n_cells += grid.second(iv).empty() ? 0 : 1;
    }
  }
  amrex::ParallelDescriptor::ReduceLongSum(n_cells);
  return n_cells;
}

amrex::Real
HaloManager::TotalMass() const
{
  amrex::Real mass = 0.0;
  for (const auto& grid : m_lines) {
    const amrex::Box& bx = grid.second.box();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      for (const auto& p : grid.second(iv)) {
        mass += p.rdata(RealData::mass);
      }
    }
  }
  amrex::ParallelDescriptor::ReduceRealSum(mass);
  return mass;
}

bool
HaloManager::Verify() const
{
  bool ok = true;
  for (const auto& grid : m_lines) {
    const amrex::Box& bx = grid.second.box();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const HaloLine& line = grid.second(iv);
      if (line.empty()) {
        continue;
      }
      ok = ok && (static_cast<int>(line.size()) == setup::n_lem);
      for (int k = 0; ok && k < static_cast<int>(line.size()); ++k) {
        const amrex::Real mass = line[k].rdata(RealData::mass);
        ok = ok && (line[k].idata(IntData::lem_index) == k) && (mass > 0.0) &&
             (std::abs(mass - line[0].rdata(RealData::mass)) <=
              1.0e-10 * line[0].rdata(RealData::mass));
      }
    }
  }
  amrex::ParallelDescriptor::ReduceBoolAnd(ok);
  return ok;
}

} // namespace clem
