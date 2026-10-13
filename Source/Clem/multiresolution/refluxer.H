/**
 * The reflux of the mass of the LEM lines at a coarse/fine interface.
 *
 * During one coarse step a coarse cell C next to the fine level exchanges
 * through its interface face twice, with different amounts: once in its own
 * splice, with the coarse transfer, and through the fine cells over their
 * substeps, with the fine transfers. The fine level is the truth, as in the
 * reflux of the LES. The difference, per interface face,
 *
 *   delta_mu_f = ( sum of the fine transfers ) - ( coarse transfer ) ,
 *
 * is a transfer still to be made, and it is made as the splice makes one,
 * between C and the coarse cell N that covers the fine cells behind the face:
 *
 *   C gave too little   delta_mu_f is cut from the HEAD of the line of C and
 *                       dropped: the fine cells already hold that mass
 *   C needs more        delta_mu_f is taken from the HEAD of the line of N
 *                       and appended at the TAIL of C. N is only copied from:
 *                       the average down that follows overwrites it
 *
 * With several interface faces, the outflux faces of C are served and the
 * arriving blocks ranked in the order of the policy of setup.H, as in the
 * splice.
 *
 * Only the mass is corrected. Nothing is done to make the species mass that C
 * exchanged match what the fine level exchanged: the composition of the mass
 * removed or added is that of the head it is taken from.
 *
 * The lines of C change in number and mass; they must be regridded after.
 *
 * Functions:
 *
 * RefluxMass : the mass reflux of the cells of a level next to the level above
 */
#ifndef CLEM_REFLUXER_H
#define CLEM_REFLUXER_H

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include "particle_container.H"

namespace clem {

// what one RefluxMass did, over the level and the ranks
struct RefluxDiagnostics
{
  amrex::Long n_cells      = 0;   // cells corrected
  amrex::Real mass_removed = 0.0; // cut from the heads and dropped
  amrex::Real mass_added   = 0.0; // copied from the covered cells
  amrex::Real max_fraction = 0.0; // max over the cells of |dm| / m_cell
};

// RefluxMass: correct the lines of the cells of level lev next to level
// lev + 1.
//   mismatch[d]  delta_mu on the faces of direction d of lev, signed along +d
//                (fine minus coarse); read on the interface faces only
//   covered      1 in the cells of lev covered by lev + 1, with one filled
//                ghost cell
//   step         step count of the level, for the policy
//   touched      on return, 1 in the cells whose line was changed
// Every line of lev must be regular on entry. Aborts if a mismatch is larger
// than the line it has to come from. On return the cell map of lev is valid.
RefluxDiagnostics RefluxMass(ClemContainer& container, int lev,
                             const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& mismatch,
                             const amrex::iMultiFab& covered, int step, amrex::iMultiFab& touched);

} // namespace clem
#endif
