#include <algorithm>

#include <AMReX_ParallelDescriptor.H>

#include "line_packing.H"

namespace clem {

// the lines of level crse under the region `need` (coarse index space),
// as a field: element k of each cell, in lem_index order, fills the
// components PackedComp(k, .). The field lives on the pieces of the coarse
// grids that `need` touches, each on the rank of its grid, so only those
// lines are packed
amrex::MultiFab
PackCoarseLines(const ClemContainer& container, const int crse, const amrex::BoxArray& need)
{
  const auto& pc                  = container.Particles();
  const amrex::BoxArray& crse_ba  = pc.ParticleBoxArray(crse);
  const auto& crse_dm             = pc.ParticleDistributionMap(crse);
  const amrex::Geometry& geom     = pc.Geom(crse);

  // `need` and its periodic images, inside the domain
  amrex::BoxList need_list;
  for (int n = 0; n < static_cast<int>(need.size()); ++n) {
    for (const auto& shift : geom.periodicity().shiftIntVect()) {
      const amrex::Box bx = amrex::shift(need[n], shift) & geom.Domain();
      if (bx.ok()) {
        need_list.push_back(bx);
      }
    }
  }
  amrex::BoxArray need_ba(need_list);
  need_ba.removeOverlap();

  // the pieces, and the coarse grid each one belongs to
  amrex::BoxList pieces;
  amrex::Vector<int> piece_rank;
  amrex::Vector<int> piece_grid;
  std::vector<std::pair<int, amrex::Box>> isects;
  for (int n = 0; n < static_cast<int>(need_ba.size()); ++n) {
    crse_ba.intersections(need_ba[n], isects);
    for (const auto& isect : isects) {
      pieces.push_back(isect.second);
      piece_rank.push_back(crse_dm[isect.first]);
      piece_grid.push_back(isect.first);
    }
  }
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pieces.isEmpty(),
                                   "FillFromCoarse: no coarse grid under the halo");

  amrex::MultiFab packed(amrex::BoxArray(pieces), amrex::DistributionMapping(piece_rank),
                         n_packed, 0);
  packed.setVal(0.0);

  std::vector<int> line;
  for (amrex::MFIter mfi(packed); mfi.isValid(); ++mfi) {
    const amrex::Box& bx  = mfi.validbox();
    const int grid_id     = piece_grid[mfi.index()];
    const auto* particles = pc.GetParticles(crse).at({grid_id, 0}).GetArrayOfStructs().data();
    const auto& out       = packed.array(mfi);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      line = container.ParticlesInCell(crse, grid_id, iv);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<int>(line.size()) == setup::n_lem,
        "FillFromCoarse: a coarse cell does not hold n_lem elements; regrid first");
      std::stable_sort(line.begin(), line.end(), [&](const int a, const int b) {
        return particles[a].idata(IntData::lem_index) < particles[b].idata(IntData::lem_index);
      });
      for (int k = 0; k < setup::n_lem; ++k) {
        for (int comp = 0; comp < RealData::ncomps; ++comp) {
          out(iv, PackedComp(k, comp)) = particles[line[k]].rdata(comp);
        }
      }
    }
  }
  return packed;
}

// the clone of the coarse line stored at coarse cell iv_crse: the same
// elements in the same order with the same state; mass and volume shared
// between the children
void
CloneLine(const amrex::Array4<const amrex::Real>& packed, const amrex::IntVect& iv_crse,
          const amrex::Real share, std::vector<ClemContainer::ParticleType>& line)
{
  line.resize(setup::n_lem);
  for (int k = 0; k < setup::n_lem; ++k) {
    auto& p = line[k];
    for (int comp = 0; comp < RealData::ncomps; ++comp) {
      p.rdata(comp) = packed(iv_crse, PackedComp(k, comp));
    }
    p.rdata(RealData::mass) *= share;
    p.rdata(RealData::vol) *= share;
    p.idata(IntData::lem_index)  = k;
    p.idata(IntData::entry_side) = 0;
  }
}


} // namespace clem
