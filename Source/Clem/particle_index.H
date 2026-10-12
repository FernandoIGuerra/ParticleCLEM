/**
 * Components carried by one CLEM particle (one element of a LEM line).
 *
 * The state is primitive: p = [rho, m, v, T, e, p, phi, Y_k], with Y_k sized
 * by the chemical mechanism (NUM_SPECIES of PelePhysics).
 */
#ifndef CLEM_PARTICLE_INDEX_H
#define CLEM_PARTICLE_INDEX_H

#include "mechanism.H"

namespace clem {

struct RealData
{
  static constexpr int rho = 0;   // density
  static constexpr int mass = 1;  // element mass
  static constexpr int vol = 2;   // element volume
  static constexpr int T = 3;     // temperature
  static constexpr int eint = 4;  // specific internal energy
  static constexpr int press = 5; // element pressure
  static constexpr int phi = 6;   // dummy passive scalar, for trivial tests
  static constexpr int Y0 = 7;    // first species mass fraction, Y_k = Y0 + k
  static constexpr int ncomps = Y0 + NUM_SPECIES;
};

struct IntData
{
  static constexpr int lem_index = 0;  // position along the 1-D LEM line
  static constexpr int entry_side = 1; // splicing entry code
  static constexpr int ncomps = 2;
};

} // namespace clem
#endif
