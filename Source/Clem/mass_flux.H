/**
 * The mass flux of the LES step, stored for the CLEM splicing.
 *
 * Each level keeps, per direction d, a face-centred field
 *
 *   mass_flux[d] = 1/2 ( F^n + F^{n+1} ) ,   F = rho* u* A   [mass / time],
 *
 * the two stages of the RK2 MOL integrator with their weight 1/2, signed
 * along +d. It is the flux of the LES density update, so per cell
 *
 *   ( rho^{n+1} - rho^n ) V = - dt sum_f +- mass_flux_f ,
 *
 * and dt mass_flux_f is the mass the lines must move through face f.
 *
 * Functions:
 *
 * AccumulateMassFlux : add one MOL stage, called from PeleC::getMOLSrcTerm
 */
#ifndef CLEM_MASS_FLUX_H
#define CLEM_MASS_FLUX_H

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_REAL.H>

namespace clem {

// add weight * flux(.., rho_comp) to mass_flux on the faces of face_box.
// flux is the extensive (area-weighted) face flux of one MOL stage, weight
// the factor of the stage in the time integrator (1/2 for RK2)
inline void
AccumulateMassFlux(const amrex::Box& face_box, const amrex::Array4<const amrex::Real>& flux,
                   const int rho_comp, const amrex::Real weight,
                   const amrex::Array4<amrex::Real>& mass_flux)
{
  amrex::LoopOnCpu(face_box, [&](int i, int j, int k) noexcept {
    mass_flux(i, j, k) += weight * flux(i, j, k, rho_comp);
  });
}

} // namespace clem
#endif
