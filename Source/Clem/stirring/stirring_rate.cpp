/**
 * stirring_rate.cpp: the rate of the stirring events of the lines.
 *
 *   StirringRate  Lambda and dt_stir of every cell of a level
 */
#include <cmath>
#include <limits>

#include <AMReX_ParallelDescriptor.H>

#include "eddy_sampling.H"
#include "setup.H"
#include "stirring_rate.H"

namespace clem {

StirringDiagnostics
StirringRate(const amrex::MultiFab& nu_t, const amrex::MultiFab& eta, amrex::MultiFab& lambda,
             amrex::MultiFab& dt_stir, const amrex::Geometry& geom, const amrex::Real dt)
{
  BL_PROFILE("clem::StirringRate()");

  constexpr amrex::Real never = std::numeric_limits<amrex::Real>::max();

  // the sizes are those of the level: the line L, the largest eddy Delta and
  // the smallest eddy the line can hold; the same range the events draw from
  const EddyRange range         = LevelEddyRange(geom);
  const amrex::Real line_length = range.line_length;
  const amrex::Real delta       = range.delta;
  const amrex::Real l_resolved  = range.l_resolved;

  const amrex::Real delta_43  = std::pow(delta, 4.0 / 3.0);
  const amrex::Real delta_m53 = std::pow(delta, -5.0 / 3.0);

  StirringDiagnostics diag;
  amrex::Real lambda_max  = 0.0;
  amrex::Real dt_stir_min = never;

  for (amrex::MFIter mfi(lambda); mfi.isValid(); ++mfi) {
    const amrex::Box& bx    = mfi.validbox();
    const auto& nu_t_arr    = nu_t.const_array(mfi);
    const auto& eta_arr     = eta.const_array(mfi);
    const auto& lambda_arr  = lambda.array(mfi);
    const auto& dt_stir_arr = dt_stir.array(mfi);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const amrex::Real eta_cell = eta_arr(iv);
      const amrex::Real l_min    = range.Smallest(eta_cell);
      const amrex::Real d_t      = nu_t_arr(iv) / setup::turbulent_schmidt;

      // no subgrid turbulence, or no inertial range left on the line
      if (!range.Stirs(eta_cell) || !(d_t > 0.0)) {
        lambda_arr(iv)  = 0.0;
        dt_stir_arr(iv) = never;
        continue;
      }

      // the third moment of f(l) over [l_min, Delta], then the rate
      const amrex::Real moment  = 1.25 * (delta_43 - std::pow(l_min, 4.0 / 3.0)) / (std::pow(l_min, -5.0 / 3.0) - delta_m53);
      const amrex::Real rate    = 27.0 * d_t / (2.0 * moment);

      lambda_arr(iv)  = rate;
      dt_stir_arr(iv) = 1.0 / (rate * line_length);

      ++diag.n_stirred;
      if (eta_cell < l_resolved) {
        ++diag.n_truncated;
      }
      lambda_max  = amrex::max(lambda_max, rate);
      dt_stir_min = amrex::min(dt_stir_min, dt_stir_arr(iv));
    }
  }

  amrex::ParallelDescriptor::ReduceLongSum(diag.n_stirred);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_truncated);
  amrex::ParallelDescriptor::ReduceRealMax(lambda_max);
  amrex::ParallelDescriptor::ReduceRealMin(dt_stir_min);
  if (diag.n_stirred > 0) {
    diag.lambda_max  = lambda_max;
    diag.dt_stir_min = dt_stir_min;
    diag.events_max  = dt / dt_stir_min;
  }
  return diag;
}

} // namespace clem
