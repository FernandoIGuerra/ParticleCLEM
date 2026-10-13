/**
 * The lines follow the grids after an AMR regrid.
 *
 * A regrid changes the grids of the levels above a base level. The LES is
 * moved to the new grids by PeleC; the lines are stored per level and per
 * grid, and are moved here. Per cell of a level whose grids changed:
 *
 *   it stays on the level    its line is refiled under the new grids and ranks
 *   it is new (refined)      it gets the clone of the line of its coarse
 *                            parent: the same elements in the same order
 *                            with the same state, mass and volume divided by
 *                            the r^d children. The children sum back to the
 *                            parent exactly
 *   it left the level        its line is dropped. The coarse cell under it
 *                            already holds the average of the lines it had
 *
 * and then every covered cell is averaged down again, so the hierarchy is as
 * after a sync. Mass, species mass and energy of the lines are conserved.
 *
 * The LES must be filled in the same way for the lines to match it: a new
 * fine cell takes the state of its parent unchanged, piecewise constant
 * (PeleC::SetNewCellsFromCoarse). A linear interpolation would give the
 * children different densities, which identical clones cannot follow.
 *
 * Requirement: every line is regular (n_lem equal-mass elements), which holds
 * between steps.
 *
 * Functions:
 *
 * AdaptToGrids : the lines of the levels above a base level onto the new grids
 */
#ifndef CLEM_GRID_ADAPTER_H
#define CLEM_GRID_ADAPTER_H

#include "particle_container.H"

namespace clem {

// what one AdaptToGrids did, over the hierarchy and the ranks
struct AdaptDiagnostics
{
  amrex::Long n_dropped = 0; // elements dropped: their cell left its level
  amrex::Long n_refined = 0; // new fine cells that got the clone of their parent
};

// the lines of the levels lbase + 1 .. finest onto their new grids. Nothing
// is done if the grids are those the lines are already filed under. On
// return every cell of every level holds n_lem elements (checked)
AdaptDiagnostics AdaptToGrids(ClemContainer& container, int lbase);

} // namespace clem
#endif
