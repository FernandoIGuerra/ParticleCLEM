/**
 * The average down of the LEM lines: the particle image of PeleC::avgDown.
 *
 * Every cell of every level holds a line. A cell covered by the finer level
 * carries the index-by-index average of the lines of its r^d children:
 * element k of the parent is the mass-weighted average of element k of each
 * child c,
 *
 *   mass_k = sum_c m_{c,k}                       (extensive: summed)
 *   vol_k  = sum_c vol_{c,k}                     rho_k = mass_k / vol_k
 *   phi_k  = sum_c m_{c,k} phi_{c,k} / mass_k    phi = eint, Y_k, phi
 *   T_k    = T( rho_k, eint_k, Y_k ) ,   press_k = p( rho_k, T_k, Y_k )
 *
 * k is the lem_index: the first element of the parent is the average of the
 * first elements of its children, and so on along the line. A child that is
 * itself covered contributes the line it got from ITS children, so a level is
 * the average of all the levels above it. Mass, volume, species mass and
 * internal energy of a covered cell are exactly those of its children, as
 * PeleC::avgDown leaves the LES. The energy, not the temperature, is the one
 * averaged: children that differ in both temperature and composition (cold
 * fuel next to hot air) do not mix to their mass-averaged temperature. The
 * temperature and the pressure are recovered from the averaged state.
 *
 * Requirement: every cell of the fine level holds a regular line, n_lem
 * elements of equal mass with lem_index 0 .. n_lem - 1 (the regridding at the
 * end of its step, and its own average down). It is checked.
 *
 * The cells of the level that are not covered are not touched.
 *
 * Functions:
 *
 * AverageDown    : the covered cells of a level from the level above
 *
 * AverageDownAll : every level, from the finest down
 */
#ifndef CLEM_AVERAGER_H
#define CLEM_AVERAGER_H

#include "particle_container.H"

namespace clem {

// what one AverageDown did, over the level and the ranks
struct AverageDiagnostics
{
  amrex::Long n_covered = 0; // covered cells of the level: lines written
  amrex::Long n_rebuilt = 0; // of them, cells whose elements had to be replaced
                             // because they did not hold a regular set of slots
};

// the lines of the covered cells of level lev from the lines of level lev + 1
// (lev < finest). On return the cell map of lev is valid
AverageDiagnostics AverageDown(ClemContainer& container, int lev);

// AverageDown for lev = finest - 1 down to 0
void AverageDownAll(ClemContainer& container);

} // namespace clem
#endif
