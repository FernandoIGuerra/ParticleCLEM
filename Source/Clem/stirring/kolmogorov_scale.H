/**
 * The Kolmogorov scale of the LES cells: the smallest eddy the stirring of
 * the lines draws.
 *
 *   eta = ( nu^3 / eps )^(1/4)
 *
 * eps is the subgrid dissipation of the LES closure, taken in equilibrium
 * with the production, eps = nu_T 2 S^d_ij S^d_ij, accumulated over the MOL
 * stages of the step from the faces (PeleC::clem_eps, PeleC::getWALELESTerm).
 * nu = mu / rho is the molecular kinematic viscosity of the cell, from its LES
 * state.
 *
 * Functions:
 *
 * KolmogorovScale : eta of every cell of a level
 */
#ifndef CLEM_KOLMOGOROV_SCALE_H
#define CLEM_KOLMOGOROV_SCALE_H

#include <AMReX_MultiFab.H>

#include "les_comps.H"
#include "particle_container.H"

namespace clem {

// what KolmogorovScale found on a level, over the cells with eps > 0
struct KolmogorovDiagnostics
{
  amrex::Long n_turbulent = 0; // cells with subgrid dissipation (eps > 0)
  amrex::Real eta_min     = 0.0;
  amrex::Real eta_max     = 0.0;
  amrex::Real delta_over_eta_max = 0.0; // the widest inertial range, Delta / eta
};

// KolmogorovScale: eta = (nu^3 / eps)^(1/4) of every cell, into `eta`, with
// nu = mu(T, Y) / rho of the LES state `les` (conserved, components `comps`)
// and eps from `eps`. A cell with eps <= 0 has no subgrid turbulence: its eta
// is 0, meaning "no inertial range, do not stir". Delta is the cell size
// along the line (x)
KolmogorovDiagnostics KolmogorovScale(const amrex::MultiFab& les, const amrex::MultiFab& eps,
                                      amrex::MultiFab& eta, const LesComps& comps,
                                      amrex::Real delta, const TransParmType* trans_parm);

} // namespace clem
#endif
