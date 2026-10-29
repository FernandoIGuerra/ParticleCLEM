/**
 * kolmogorov_scale.cpp: the Kolmogorov scale of the LES cells.
 *
 *   KolmogorovScale  eta = (nu^3 / eps)^(1/4) of every cell of a level
 */
#include <cmath>
#include <limits>

#include <AMReX_ParallelDescriptor.H>

#include "PelePhysics.H"

#include "kolmogorov_scale.H"

namespace clem {

KolmogorovDiagnostics
KolmogorovScale(const amrex::MultiFab& les, const amrex::MultiFab& eps, amrex::MultiFab& eta,
                const LesComps& comps, const amrex::Real delta, const TransParmType* trans_parm)
{
  BL_PROFILE("clem::KolmogorovScale()");

  auto trans = pele::physics::PhysicsType::transport();

  // only the viscosity
  constexpr bool get_xi    = false;
  constexpr bool get_mu    = true;
  constexpr bool get_lam   = false;
  constexpr bool get_ddiag = false;
  constexpr bool get_chi   = false;

  KolmogorovDiagnostics diag;
  amrex::Real eta_min = std::numeric_limits<amrex::Real>::max();
  amrex::Real eta_max = 0.0;

  for (amrex::MFIter mfi(eta); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const auto& state    = les.const_array(mfi);
    const auto& eps_arr  = eps.const_array(mfi);
    const auto& eta_arr  = eta.array(mfi);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const amrex::Real eps_cell = eps_arr(iv);
      if (!(eps_cell > 0.0)) {
        eta_arr(iv) = 0.0;
        continue;
      }

      const amrex::Real rho  = state(iv, comps.rho);
      amrex::Real temp       = state(iv, comps.temp);
      amrex::Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = state(iv, comps.spec + n) / rho;
      }

      amrex::Real d_diag[NUM_SPECIES];
      amrex::Real mu       = 0.0;
      amrex::Real xi       = 0.0;
      amrex::Real lam      = 0.0;
      amrex::Real* chi_mix = nullptr;
      trans.transport(get_xi, get_mu, get_lam, get_ddiag, get_chi, temp, rho, Y, d_diag, chi_mix,
                      mu, xi, lam, trans_parm);

      const amrex::Real nu       = mu / rho;
      const amrex::Real eta_cell = std::pow(nu * nu * nu / eps_cell, 0.25);
      eta_arr(iv)                = eta_cell;

      ++diag.n_turbulent;
      eta_min = amrex::min(eta_min, eta_cell);
      eta_max = amrex::max(eta_max, eta_cell);
    }
  }

  amrex::ParallelDescriptor::ReduceLongSum(diag.n_turbulent);
  amrex::ParallelDescriptor::ReduceRealMin(eta_min);
  amrex::ParallelDescriptor::ReduceRealMax(eta_max);
  if (diag.n_turbulent > 0) {
    diag.eta_min            = eta_min;
    diag.eta_max            = eta_max;
    diag.delta_over_eta_max = delta / eta_min;
  }
  return diag;
}

} // namespace clem
