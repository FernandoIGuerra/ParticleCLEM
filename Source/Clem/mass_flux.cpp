/**
 * mass_flux.cpp: the check of the stored mass flux against the LES.
 *
 *   PeleC::CheckClemMassFlux  per level, the largest relative defect of
 *                             ( rho^{n+1} - rho^n ) V + dt div( mass_flux )
 */
#include <AMReX_ParallelDescriptor.H>

#include "PeleC.H"
#include "IndexDefines.H"

void
PeleC::CheckClemMassFlux(const amrex::Real dt)
{
  BL_PROFILE("PeleC::CheckClemMassFlux()");

  const amrex::MultiFab& state_old = get_old_data(State_Type);
  const amrex::MultiFab& state_new = get_new_data(State_Type);
  const auto cell_size             = geom.CellSizeArray();
  const amrex::Real cell_vol       = AMREX_D_TERM(cell_size[0], *cell_size[1], *cell_size[2]);

  amrex::Real max_defect = 0.0;
  for (amrex::MFIter mfi(state_new); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const auto& rho_old  = state_old.const_array(mfi, URHO);
    const auto& rho_new  = state_new.const_array(mfi, URHO);
    AMREX_D_TERM(const auto& flux_x = clem_mass_flux[0].const_array(mfi);
                 , const auto& flux_y = clem_mass_flux[1].const_array(mfi);
                 , const auto& flux_z = clem_mass_flux[2].const_array(mfi);)

    amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
      // net mass rate leaving the cell
      const amrex::Real net_out = AMREX_D_TERM(flux_x(i + 1, j, k) - flux_x(i, j, k),
                                               +flux_y(i, j + 1, k) - flux_y(i, j, k),
                                               +flux_z(i, j, k + 1) - flux_z(i, j, k));
      const amrex::Real defect = (rho_new(i, j, k) - rho_old(i, j, k)) * cell_vol + dt * net_out;
      max_defect = amrex::max(max_defect, std::abs(defect) / (rho_old(i, j, k) * cell_vol));
    });
  }
  amrex::ParallelDescriptor::ReduceRealMax(max_defect);

  amrex::Print() << "CLEM: level " << level << " mass flux defect, max |d(rho V) + dt div F| / (rho V) = "
                 << max_defect << std::endl;
}
