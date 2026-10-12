/**
 * Small helpers on the equation of state, shared by the operators of the lines.
 *
 * Functions:
 *
 * TemperatureFromEnergy : T from (rho, e, Y)
 */
#ifndef CLEM_EOS_UTIL_H
#define CLEM_EOS_UTIL_H

#include <cmath>

#include <AMReX_REAL.H>

#include "PelePhysics.H"

namespace clem {

// T from (rho, e, Y), starting from the guess in temp. The inversion of the
// EOS is polished with a few Newton steps on e(T), so that e(rho, T, Y) is
// the given energy to round-off
inline void
TemperatureFromEnergy(const amrex::Real rho, const amrex::Real eint, const amrex::Real* Y,
                      amrex::Real& temp)
{
  auto eos = pele::physics::PhysicsType::eos();
  eos.REY2T(rho, eint, Y, temp);
  for (int iter = 0; iter < 3; ++iter) {
    amrex::Real eint_T = 0.0;
    amrex::Real cv     = 0.0;
    eos.RTY2E(rho, temp, Y, eint_T);
    eos.RTY2Cv(rho, temp, Y, cv);
    const amrex::Real dT = (eint - eint_T) / cv;
    temp += dT;
    if (std::abs(dT) < 1.0e-12 * temp) {
      break;
    }
  }
}

} // namespace clem
#endif
