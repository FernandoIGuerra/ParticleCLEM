/**
 * This is a derived class from Regridder.
 *
 * EqualMassRegridder: every line back to setup::n_lem elements of EQUAL MASS
 * m_c / n_lem, m_c the mass of the line (Maxwell thesis, p. 172).
 *
 * The line is walked in lem_index order and its mass poured, in that order,
 * into the n_lem new elements. A new element j that received the pieces dm_i
 * of the old elements i gets
 *
 *   mass_j = sum_i dm_i = m_c / n_lem
 *   vol_j  = sum_i dm_i / rho_i               rho_j = mass_j / vol_j
 *   phi_j  = sum_i dm_i phi_i / mass_j        phi = eint, Y_k, phi
 *   T_j    = T( rho_j, eint_j, Y_j ) ,   press_j = p( rho_j, T_j, Y_j )
 *
 * so the mass, the volume, the species mass and the internal energy of the
 * line are conserved. The energy, not the temperature, is the one averaged:
 * pieces that differ in both temperature and composition do not mix to their
 * mass-averaged temperature. The temperature and the pressure are recovered
 * from the merged state, so the element stays consistent with the equation
 * of state, which the constant-pressure diffusion and the isentropic
 * pressure change need.
 *
 * Steps:
 *   1. a cell with fewer than n_lem elements gets empty ones at the end of
 *      its line, so it has a slot for every new element
 *   2. the merge, per cell, into the first n_lem slots; the other slots are
 *      invalidated
 *   3. Redistribute of the level removes the invalidated elements
 */
#ifndef CLEM_EQUAL_MASS_REGRIDDER_H
#define CLEM_EQUAL_MASS_REGRIDDER_H

#include "regridder.H"

namespace clem {

class EqualMassRegridder : public Regridder<EqualMassRegridder>
{
  friend Regridder;

protected:
  RegridDiagnostics RegridImplementation(ClemContainer& container, int lev) const;

private:
  // step 1; returns the number of elements created on this rank
  amrex::Long TopUpLines(ClemContainer& container, int lev) const;

  // step 2; returns the number of elements invalidated on this rank
  amrex::Long MergeLines(ClemContainer& container, int lev) const;

public:
  static constexpr const char* name = "equal mass, e and Y mass-averaged, T and p from the EOS";

  EqualMassRegridder()                                     = default;
  ~EqualMassRegridder()                                    = default;
  EqualMassRegridder(EqualMassRegridder const&)            = delete;
  EqualMassRegridder& operator=(EqualMassRegridder const&) = delete;
  EqualMassRegridder(EqualMassRegridder&&)                 = delete;
  EqualMassRegridder& operator=(EqualMassRegridder&&)      = delete;
};

} // namespace clem
#endif
