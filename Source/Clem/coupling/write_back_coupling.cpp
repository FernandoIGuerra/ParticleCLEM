#include <cmath>

#include <AMReX_ParallelDescriptor.H>

#include "PelePhysics.H"

#include "write_back_coupling.H"

namespace clem {

namespace {
// one element compressed isentropically, at frozen composition, to the
// density rho_new and temperature temp_new; its pressure becomes press_new.
// The mass is fixed, so the volume follows; the energy changes by the work
template <typename ParticleType>
void
CompressElement(ParticleType& p, const amrex::Real* Y, const amrex::Real rho_new,
                const amrex::Real temp_new, const amrex::Real press_new)
{
  auto eos = pele::physics::PhysicsType::eos();

  // difference form: e += e(new) - e(old)
  amrex::Real eint_old = 0.0;
  amrex::Real eint_new = 0.0;
  eos.RTY2E(p.rdata(RealData::rho), p.rdata(RealData::T), Y, eint_old);
  eos.RTY2E(rho_new, temp_new, Y, eint_new);

  p.rdata(RealData::rho)   = rho_new;
  p.rdata(RealData::vol)   = p.rdata(RealData::mass) / rho_new;
  p.rdata(RealData::T)     = temp_new;
  p.rdata(RealData::press) = press_new;
  p.rdata(RealData::eint) += eint_new - eint_old;
}

amrex::Real
CellVolume(const ClemContainer& container, const int lev)
{
  const auto cell_size = container.Particles().Geom(lev).CellSizeArray();
  return AMREX_D_TERM(cell_size[0], *cell_size[1], *cell_size[2]);
}
} // namespace

void
ConservedSnapshot(const ClemContainer& container, const int lev, amrex::MultiFab& conserved)
{
  BL_PROFILE("clem::ConservedSnapshot()");
  const auto& pc = container.Particles();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    conserved.boxArray() == pc.ParticleBoxArray(lev) &&
      conserved.DistributionMap() == pc.ParticleDistributionMap(lev),
    "ConservedSnapshot: conserved must live on the particle grids of this level");
  AMREX_ALWAYS_ASSERT(conserved.nComp() >= NUM_SPECIES + 1);

  const amrex::Real inv_vol = 1.0 / CellVolume(container, lev);

  // straight from the particles, so the result does not depend on the cell map
  conserved.setVal(0.0);
  for (ClemParConstIter pti(pc, lev); pti.isValid(); ++pti) {
    const auto& aos = pti.GetArrayOfStructs();
    const auto& sum = conserved.array(pti.index());

    for (int n = 0; n < pti.numParticles(); ++n) {

      const auto& p             = aos[n];
      const amrex::IntVect iv   = pc.Index(p, lev);
      const amrex::Real weight  = p.rdata(RealData::mass) * inv_vol;

      for (int k = 0; k < NUM_SPECIES; ++k) {
        sum(iv, k) += weight * p.rdata(RealData::Y0 + k);
      }
      sum(iv, NUM_SPECIES) += weight * p.rdata(RealData::eint);
    }
  }
}

void
WriteBackCoupling::ForwardCouplingImplementation(ClemContainer& container, const int lev,
                                                 const amrex::MultiFab& les_new,
                                                 const LesComps& comps) const
{
  BL_PROFILE("clem::WriteBackCoupling::ForwardCoupling()");
  auto& pc = container.Particles();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    les_new.boxArray() == pc.ParticleBoxArray(lev) &&
      les_new.DistributionMap() == pc.ParticleDistributionMap(lev),
    "ForwardCoupling: LES state must live on the particle grids of this level");

  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev);
  }
  auto eos = pele::physics::PhysicsType::eos();

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    const auto& state    = les_new.const_array(grid_id);
    auto* particles      = pti.GetArrayOfStructs().data();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      // the target: the pressure of the LES cell, once per cell
      const amrex::Real rho_les = state(iv, comps.rho);

      amrex::Real Y_les[NUM_SPECIES];
      for (int k = 0; k < NUM_SPECIES; ++k) {
        Y_les[k] = state(iv, comps.spec + k) / rho_les;
      }

      amrex::Real press_les = 0.0;
      eos.RTY2P(rho_les, state(iv, comps.temp), Y_les, press_les);

      for (const int n : container.ParticlesInCell(lev, grid_id, iv)) {
        auto& p                = particles[n];
        const amrex::Real rho  = p.rdata(RealData::rho);
        const amrex::Real temp = p.rdata(RealData::T);

        if (rho <= 0.0 || p.rdata(RealData::press) <= 0.0) {
          continue;
        }

        const amrex::Real* Y = &p.rdata(RealData::Y0);

        // the local isentropic exponent of the element
        amrex::Real gamma = 0.0;
        eos.RTY2G(rho, temp, Y, gamma);

        const amrex::Real press_ratio = press_les / p.rdata(RealData::press);
        CompressElement(p, Y, rho * std::pow(press_ratio, 1.0 / gamma),
                        temp * std::pow(press_ratio, (gamma - 1.0) / gamma), press_les);
      }
    }
  }
}

amrex::Real
WriteBackCoupling::RenormalizeVolume(ClemContainer& container, const int lev) const
{
  auto& pc = container.Particles();
  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev);
  }
  auto eos                   = pele::physics::PhysicsType::eos();
  const amrex::Real cell_vol = CellVolume(container, lev);

  amrex::Real max_defect = 0.0;
  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    auto* particles      = pti.GetArrayOfStructs().data();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& in_cell = container.ParticlesInCell(lev, grid_id, iv);

      amrex::Real line_vol = 0.0;
      for (const int n : in_cell) {
        line_vol += particles[n].rdata(RealData::vol);
      }
      if (line_vol <= 0.0) {
        continue;
      }
      // the same ratio for every element of the cell
      const amrex::Real f = cell_vol / line_vol;
      max_defect          = amrex::max<amrex::Real>(max_defect, std::abs(f - 1.0));
      if (f == 1.0) {
        continue;
      }

      for (const int n : in_cell) {
        auto& p                = particles[n];
        const amrex::Real rho  = p.rdata(RealData::rho);
        const amrex::Real temp = p.rdata(RealData::T);
        if (rho <= 0.0 || p.rdata(RealData::press) <= 0.0) {
          continue;
        }
        const amrex::Real* Y = &p.rdata(RealData::Y0);

        amrex::Real gamma = 0.0;
        eos.RTY2G(rho, temp, Y, gamma);

        CompressElement(p, Y, rho / f, temp * std::pow(f, -(gamma - 1.0)),
                        p.rdata(RealData::press) * std::pow(f, -gamma));
      }
    }
  }
  return max_defect;
}

CouplingDiagnostics
WriteBackCoupling::BackwardCouplingImplementation(ClemContainer& container, const int lev,
                                                  amrex::MultiFab& les_new,
                                                  const amrex::MultiFab& les_old,
                                                  const amrex::MultiFab& conserved_pre,
                                                  const LesComps& comps) const
{
  BL_PROFILE("clem::WriteBackCoupling::BackwardCoupling()");
  AMREX_ALWAYS_ASSERT(conserved_pre.nComp() >= NUM_SPECIES + 1);

  CouplingDiagnostics diag;

  // 1. the lines back into the volume of their cells
  diag.max_volume_defect = RenormalizeVolume(container, lev);

  // 2. the conserved filter after everything the lines did in the step
  amrex::MultiFab conserved_post(les_new.boxArray(), les_new.DistributionMap(), NUM_SPECIES + 1,
                                 0);
  ConservedSnapshot(container, lev, conserved_post);

  for (amrex::MFIter mfi(les_new); mfi.isValid(); ++mfi) {
    const amrex::Box& bx  = mfi.validbox();
    const auto& state_new = les_new.array(mfi);
    const auto& state_old = les_old.const_array(mfi);
    const auto& pre       = conserved_pre.const_array(mfi);
    const auto& post      = conserved_post.const_array(mfi);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const amrex::Real rho_new = state_new(iv, comps.rho);

      // 3. ( rho Y_k )^n + Delta[rho Y_k], floored, then projected onto the
      // density of the LES
      amrex::Real rho_Y[NUM_SPECIES];
      amrex::Real rho_Y_sum = 0.0;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        rho_Y[k] = amrex::max<amrex::Real>( 0.0, state_old(iv, comps.spec + k) + (post(iv, k) - pre(iv, k)));
        rho_Y_sum += rho_Y[k];
      }
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rho_new > 0.0 && rho_Y_sum > 0.0, "BackwardCoupling: a cell holds no mass");
      
      const amrex::Real projection = rho_new / rho_Y_sum;
      diag.max_projection = amrex::max<amrex::Real>(diag.max_projection, std::abs(projection - 1.0));
      
      for (int k = 0; k < NUM_SPECIES; ++k) {
        state_new(iv, comps.spec + k) = rho_Y[k] * projection;
      }

      // 4. ( rho e )^n + Delta[rho e]; rho E rebuilt with the kinetic energy
      // of the LES, whose momentum is not touched
      amrex::Real kinetic_old = 0.0;
      amrex::Real kinetic_new = 0.0;

      for (int d = 0; d < 3; ++d) {
        kinetic_old += state_old(iv, comps.mom + d) * state_old(iv, comps.mom + d);
        kinetic_new += state_new(iv, comps.mom + d) * state_new(iv, comps.mom + d);
      }
      
      kinetic_old *= 0.5 / state_old(iv, comps.rho);
      kinetic_new *= 0.5 / rho_new;

      const amrex::Real rho_e = (state_old(iv, comps.eden) - kinetic_old) +
                                (post(iv, NUM_SPECIES) - pre(iv, NUM_SPECIES));
      state_new(iv, comps.eint) = rho_e;
      state_new(iv, comps.eden) = rho_e + kinetic_new;
    }
  }

  amrex::ParallelDescriptor::ReduceRealMax(diag.max_volume_defect);
  amrex::ParallelDescriptor::ReduceRealMax(diag.max_projection);
  return diag;
}

} // namespace clem
