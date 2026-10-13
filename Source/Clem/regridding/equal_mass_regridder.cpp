#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_ParallelDescriptor.H>

#include "PelePhysics.H"

#include "equal_mass_regridder.H"
#include "eos_util.H"
#include "setup.H"

namespace clem {

namespace {
using ParticleType = ClemContainer::ParticleType;

} // namespace

RegridDiagnostics
EqualMassRegridder::RegridImplementation(ClemContainer& container, const int lev,
                                         const amrex::iMultiFab* only_cells) const
{
  BL_PROFILE("clem::EqualMassRegridder::Regrid()");
  AMREX_ALWAYS_ASSERT(container.IsDefined());

  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev);
  }

  RegridDiagnostics diag;

  // 1. the empty elements go straight into the tile of their cell, so only
  // the map has to see them
  diag.n_added = TopUpLines(container, lev, only_cells);

  amrex::ParallelDescriptor::ReduceLongSum(diag.n_added);
  if (diag.n_added > 0) {
    container.RebuildCellMap(lev);
  }

  // 2.
  diag.n_removed = MergeLines(container, lev, only_cells);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_removed);

  // 3. only if some rank invalidated elements; otherwise the storage, and
  // with it the map, did not change
  if (diag.n_removed > 0) {
    container.Redistribute(lev);
    container.RebuildCellMap(lev);
  }
  return diag;
}

amrex::Long
EqualMassRegridder::TopUpLines(ClemContainer& container, const int lev,
                               const amrex::iMultiFab* only_cells) const
{
  auto& pc             = container.Particles();
  const auto prob_lo   = pc.Geom(lev).ProbLoArray();
  const auto cell_size = pc.Geom(lev).CellSizeArray();

  amrex::Long n_added = 0;
  for (amrex::MFIter mfi = pc.MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const int grid_id    = mfi.index();
    auto& tile           = pc.DefineAndReturnParticleTile(lev, grid_id, mfi.LocalTileIndex());

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (only_cells != nullptr && only_cells->const_array(mfi)(iv) != 1) {
        continue;
      }
      // only the count is read: push_back may move the storage of the tile
      const int n_in_cell = static_cast<int>(container.ParticlesInCell(lev, grid_id, iv).size());

      for (int k = n_in_cell; k < setup::n_lem; ++k) {
        ParticleType p;
        p.id()  = ParticleType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();

        for (int comp = 0; comp < RealData::ncomps; ++comp) {
          p.rdata(comp) = 0.0;
        }
        // at the end of the line: the merge fills the slots in line order
        PlaceOnLine(p, iv, k, setup::n_lem, prob_lo, cell_size);
        p.idata(IntData::lem_index)  = std::numeric_limits<int>::max() - setup::n_lem + k;
        p.idata(IntData::entry_side) = 0;
        tile.push_back(p);
        ++n_added;
      }
    }
  }
  return n_added;
}

amrex::Long
EqualMassRegridder::MergeLines(ClemContainer& container, const int lev,
                               const amrex::iMultiFab* only_cells) const
{
  auto& pc             = container.Particles();
  const auto prob_lo   = pc.Geom(lev).ProbLoArray();
  const auto cell_size = pc.Geom(lev).CellSizeArray();
  auto eos             = pele::physics::PhysicsType::eos();

  constexpr int ncomps = RealData::ncomps;
  // reused for every cell: the old elements of the line, and its slots in
  // line order
  std::vector<amrex::Real> old_data;
  std::vector<int> slots;

  amrex::Long n_removed = 0;
  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    auto* particles      = pti.GetArrayOfStructs().data();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (only_cells != nullptr && only_cells->const_array(grid_id)(iv) != 1) {
        continue;
      }
      const auto& in_cell = container.ParticlesInCell(lev, grid_id, iv);

      slots.assign(in_cell.begin(), in_cell.end());
      
      std::sort(slots.begin(), slots.end(), [&](const int a, const int b) {
        return particles[a].idata(IntData::lem_index) < particles[b].idata(IntData::lem_index);
      });

      const int n_old = static_cast<int>(slots.size());
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_old >= setup::n_lem, "CLEM regrid: a line has fewer slots than n_lem");

      // the old elements; their slots are reused as the new ones
      old_data.resize(n_old * ncomps);
      amrex::Real line_mass = 0.0;

      for (int i = 0; i < n_old; ++i) {
        const auto& p = particles[slots[i]];
        for (int comp = 0; comp < ncomps; ++comp) {
          old_data[i * ncomps + comp] = p.rdata(comp);
        }
        line_mass += p.rdata(RealData::mass);
      }

      // the slots beyond n_lem leave the line
      for (int i = setup::n_lem; i < n_old; ++i) {
        particles[slots[i]].id() = -particles[slots[i]].id();
        ++n_removed;
      }

      for (int j = 0; j < setup::n_lem; ++j) {
        auto& p = particles[slots[j]];
        PlaceOnLine(p, iv, j, setup::n_lem, prob_lo, cell_size);
        p.idata(IntData::entry_side) = 0;
      }
      if (line_mass <= 0.0) {
        continue; // an empty line: n_lem empty elements
      }

      const amrex::Real new_mass = line_mass / static_cast<amrex::Real>(setup::n_lem);

      for (int j = 0; j < setup::n_lem; ++j) {
        auto& p = particles[slots[j]];
        for (int comp = 0; comp < ncomps; ++comp) {
          p.rdata(comp) = 0.0;
        }
      }

      // pour the old elements, in line order, into the new ones: the new
      // element j is closed when it holds new_mass, the old element i is
      // left when it is empty. e, phi and Y_k accumulate as dm phi, the
      // volume as dm / rho of the old element
      int i = 0;
      int j = 0;
      while (i < n_old && j < setup::n_lem) {
        auto& p                 = particles[slots[j]];
        amrex::Real& old_mass   = old_data[i * ncomps + RealData::mass];
        const amrex::Real need  = new_mass - p.rdata(RealData::mass);
        // the last new element takes whatever is left, so round-off cannot
        // leave mass behind
        const bool close      = (j < setup::n_lem - 1) && (old_mass >= need);
        const amrex::Real dm  = close ? need : old_mass;

        if (dm > 0.0) {
          p.rdata(RealData::mass)   += dm;
          p.rdata(RealData::vol)    += dm / old_data[i * ncomps + RealData::rho];
          p.rdata(RealData::eint)   += dm * old_data[i * ncomps + RealData::eint];
          p.rdata(RealData::T)      += dm * old_data[i * ncomps + RealData::T];
          p.rdata(RealData::phi)    += dm * old_data[i * ncomps + RealData::phi];
          
          for (int k = 0; k < NUM_SPECIES; ++k) {
            p.rdata(RealData::Y0 + k) += dm * old_data[i * ncomps + RealData::Y0 + k];
          }
          old_mass -= dm;
        }
        if (close) {
          ++j;
        } else {
          ++i;
        }
      }

      // mass averages, rho from the conserved mass and volume, and the
      // temperature and the pressure from the equation of state
      for (j = 0; j < setup::n_lem; ++j) {
        auto& p                = particles[slots[j]];
        const amrex::Real mass = p.rdata(RealData::mass);
        
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          mass > 0.0 && p.rdata(RealData::vol) > 0.0 &&std::abs(mass - new_mass) <= 1.0e-10 * new_mass,
          "CLEM regrid: a new element did not get its equal share of the mass");

        const amrex::Real inv_mass  = 1.0 / mass;
        p.rdata(RealData::eint)     *= inv_mass;
        p.rdata(RealData::T)        *= inv_mass; // only the guess of the inversion below
        p.rdata(RealData::phi)      *= inv_mass;

        for (int k = 0; k < NUM_SPECIES; ++k) {
          p.rdata(RealData::Y0 + k) *= inv_mass;
        }

        p.rdata(RealData::rho)      = mass / p.rdata(RealData::vol);

        const amrex::Real* Y = &p.rdata(RealData::Y0);
        
        // the energy is the merged one; the temperature follows from it
        TemperatureFromEnergy(p.rdata(RealData::rho), p.rdata(RealData::eint), Y, p.rdata(RealData::T));
        eos.RTY2P(p.rdata(RealData::rho), p.rdata(RealData::T), Y, p.rdata(RealData::press));
      }
    }
  }
  return n_removed;
}

} // namespace clem
