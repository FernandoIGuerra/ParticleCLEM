/**
 * Where the LES state keeps its fields.
 */
#ifndef CLEM_LES_COMPS_H
#define CLEM_LES_COMPS_H

namespace clem {

// rho, rho u (first of 3), rho E, rho e, T, rho Y_k (first of NUM_SPECIES)
struct LesComps
{
  int rho;
  int mom;
  int eden;
  int eint;
  int temp;
  int spec;
};

} // namespace clem
#endif
