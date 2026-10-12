/**
 * The halo (ghost) lines of a fine level at a coarse/fine interface.
 *
 * A fine cell next to the interface has no fine neighbour across it: the
 * fluid there is carried by the line of a COARSE cell. To give the fine
 * level a neighbour of its own size, every fine cell position of the first
 * layer around the fine grids that lies on the coarser level gets a halo
 * line: an exact copy of the line of the coarse cell C that contains it,
 *
 *   halo line  =  the n_lem elements of C, in the same order (lem_index),
 *                 with the same state; only mass and vol are divided by the
 *                 r^d children of C (4 in 2D, 8 in 3D for r = 2)
 *
 * so the r^d children of C sum back to the line of C exactly. The coarse
 * lines must be regridded: every cell holds n_lem elements of equal mass.
 *
 * The halo lines are NOT particles of the container (a particle placed in a
 * halo cell would be filed on the coarse level): they are a copy, held here
 * per fine grid, refilled whenever the coarse lines change. Nothing in the
 * container is modified.
 *
 * A halo position covered by another fine grid has a real fine line, and one
 * outside a non-periodic domain is a boundary: neither gets a halo line.
 *
 * HaloManager holds the halo of ONE fine level.
 *
 * Functions:
 *
 * FillFromCoarse : build the halo lines of a level from the level below
 *
 * Line           : the halo line of one cell
 *
 * Verify         : every halo line holds n_lem equal-mass elements, in order
 */
#ifndef CLEM_HALO_MANAGER_H
#define CLEM_HALO_MANAGER_H

#include <map>
#include <vector>

#include <AMReX_BaseFab.H>
#include <AMReX_IntVect.H>

#include "particle_container.H"

namespace clem {

class HaloManager
{
public:
  // the line of one halo cell, head first; empty where there is no halo. The
  // elements are particles so that the operators of the lines apply to them;
  // their position and id are not used
  using HaloLine = std::vector<ClemContainer::ParticleType>;

  HaloManager() = default;

  // FillFromCoarse: build the halo lines of level lev (lev >= 1) from the
  // lines of level lev - 1, on one layer of fine cells around every grid:
  //   1. the coarse lines under that layer are packed, element by element in
  //      lem_index order, into a field on the coarse grids
  //   2. that field is copied onto the coarsened layer of each fine grid: this
  //      is what crosses ranks and periodic boundaries
  //   3. each halo cell clones the line of its coarse parent, mass and volume
  //      divided by r^d
  // The cell map of lev - 1 is rebuilt if stale. Any previous halo is dropped.
  void FillFromCoarse(ClemContainer& container, int lev);

  void Clear();
  bool IsDefined() const { return m_lev >= 0; }
  int Level() const { return m_lev; }

  // is fine cell iv a halo cell of grid grid_id (a grid of this rank)
  bool IsHaloCell(int grid_id, const amrex::IntVect& iv) const;

  // the halo line of cell iv of grid grid_id; empty if iv is not a halo cell
  const HaloLine& Line(int grid_id, const amrex::IntVect& iv) const;

  // number of halo cells and their mass, all ranks
  amrex::Long NumHaloCells() const;
  amrex::Real TotalMass() const;

  // as filled: every halo line holds n_lem elements of equal mass, in
  // lem_index order 0 .. n_lem - 1 (all ranks)
  bool Verify() const;

private:
  int m_lev = -1;
  // per fine grid of this rank, on its box grown by one cell
  std::map<int, amrex::BaseFab<HaloLine>> m_lines;
};

} // namespace clem
#endif
