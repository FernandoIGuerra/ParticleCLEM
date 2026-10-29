/**
 * The rate of the stirring events of the lines (Chakravarthy & Menon 2001).
 *
 * The triplet maps of a line draw their size l from the inertial range,
 *
 *   f(l) = (5/3) l^(-8/3) / ( l_min^(-5/3) - Delta^(-5/3) ),   l_min <= l <= Delta
 *
 * and their rate per unit length Lambda is set so that the stirring
 * reproduces the turbulent diffusivity of the LES closure,
 *
 *   D_T    = nu_T / Sc_T = (2/27) Lambda I
 *   I      = Int l^3 f(l) dl = (5/4) (Delta^(4/3) - l_min^(4/3)) / (l_min^(-5/3) - Delta^(-5/3))
 *   Lambda = 27 D_T / (2 I)                                     [1/(cm s)]
 *
 * A line of length L sees Lambda L events per unit time, one every
 *
 *   dt_stir = 1 / (Lambda L)                                    [s]
 *
 * on average. The sizes:
 *
 * (the sizes come from LevelEddyRange, eddy_sampling.H, shared with the draws
 * of the events)
 *
 *   L      the line, the cell size along x
 *   Delta  the largest eddy, min(filter size, L); the filter size is
 *          (dx dy dz)^(1/3), as in the WALE closure
 *   l_min  the smallest eddy, max(eta, n_min l_particle), l_particle = L / n_lem
 *          the length of an element: a triplet map needs n_min =
 *          setup::stirring_min_elements elements, so the line cannot hold
 *          eddies below that, and the range is truncated there; Lambda is
 *          computed over the truncated range, so the eddies the line holds
 *          still give the D_T of the LES
 *
 * A cell with no inertial range, eta = 0 (eps = 0, no subgrid turbulence) or
 * l_min >= Delta, is not stirred: Lambda = 0 and dt_stir = the largest Real,
 * an event that never comes.
 *
 * Functions:
 *
 * StirringRate : Lambda and dt_stir of every cell of a level
 */
#ifndef CLEM_STIRRING_RATE_H
#define CLEM_STIRRING_RATE_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace clem {

// what StirringRate found on a level, over the stirred cells
struct StirringDiagnostics
{
  amrex::Long n_stirred   = 0;   // cells with an inertial range (Lambda > 0)
  amrex::Long n_truncated = 0;   // of them, those with eta below the line resolution
  amrex::Real lambda_max  = 0.0; // the largest event rate per unit length [1/(cm s)]
  amrex::Real dt_stir_min = 0.0; // the shortest mean time between events of a line [s]
  amrex::Real events_max  = 0.0; // the most events of a line in this LES step, dt / dt_stir
};

// StirringRate: from nu_T (`nu_t`) and the Kolmogorov scale (`eta`, 0 for a
// cell with no subgrid turbulence) of every cell of a level, write the event
// rate per unit length into `lambda` and the mean time between the events of
// its line into `dt_stir`. geom gives the cell sizes, dt the LES step (for
// the diagnostics only)
StirringDiagnostics StirringRate(const amrex::MultiFab& nu_t, const amrex::MultiFab& eta,
                                 amrex::MultiFab& lambda, amrex::MultiFab& dt_stir,
                                 const amrex::Geometry& geom, amrex::Real dt);

} // namespace clem
#endif
