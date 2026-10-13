/**
 * The lines of a level as a field, and the clone of a line onto its children.
 *
 * Shared by the operators that copy coarse lines onto fine cells: the halo
 * (halo_manager.H) and the lines of the cells that an AMR regrid has just
 * refined (grid_adapter.H).
 *
 * A regular line (n_lem elements) is packed element by element, in lem_index
 * order: element k fills the components PackedComp(k, .), one per real datum
 * of the particle.
 *
 * Functions:
 *
 * PackCoarseLines : the lines of a level under a region, as a field
 *
 * CloneLine       : the line of a coarse cell, with a share of its mass
 */
#ifndef CLEM_LINE_PACKING_H
#define CLEM_LINE_PACKING_H

#include <vector>

#include <AMReX_MultiFab.H>

#include "particle_container.H"
#include "setup.H"

namespace clem {

// components of the packed field: element k of the line, real datum comp
constexpr int n_packed = setup::n_lem * RealData::ncomps;
inline int PackedComp(const int k, const int comp) { return k * RealData::ncomps + comp; }

// the lines of level crse under the region `need` (coarse index space). The
// field lives on the pieces of the coarse grids that `need` touches, each on
// the rank of its grid, so only those lines are packed. Every cell packed
// must hold n_lem elements. The cell map of crse must be valid
amrex::MultiFab PackCoarseLines(const ClemContainer& container, int crse,
                                const amrex::BoxArray& need);

// the clone of the coarse line stored at coarse cell iv_crse of `packed`: the
// same elements in the same order with the same state; mass and volume
// multiplied by `share`
void CloneLine(const amrex::Array4<const amrex::Real>& packed, const amrex::IntVect& iv_crse,
               amrex::Real share, std::vector<ClemContainer::ParticleType>& line);

} // namespace clem
#endif
