/**
 * The CLEM particle container.
 *
 * ClemParticleContainer is the AMReX storage: one AmrParticleContainer over
 * the whole hierarchy, every component inside the particle (AoS).
 *
 * ClemContainer wraps it: it owns the storage and is the place where the
 * operations on the LEM lines are declared.
 *
 * The lines live on EVERY level: each cell of each level holds its own line,
 * covered or not, so the particles are block-structured as the LES mesh.
 *
 * The operations work on ONE level (argument lev); the form without lev loops
 * over the levels.
 *
 * Functions:
 *
 * InitializeFromLES  : seed the lines of a level from the LES state
 *
 * GetFilteredScalars : the filtered (cell) values of the lines of a level
 *
 * ParticlesInCell    : the particles of one cell, through the cell map
 *
 * Redistribute       : AMReX Redistribute of a level, its cell map becomes stale
 *
 * DiffusionStirringReactor : the lines of a level advanced over one LES step
 *
 * WritePlotFile      : every particle of every level
 *
 * WriteLeafPlotFile  : only the particles of the cells not covered by a finer level
 */
#ifndef CLEM_PARTICLE_CONTAINER_H
#define CLEM_PARTICLE_CONTAINER_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <AMReX_AmrCore.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_BaseFab.H>
#include <AMReX_MultiFab.H>

#include "PelePhysics.H"

#include "particle_index.H"

namespace clem {

// the transport coefficients of PelePhysics; PeleC owns them and passes them in
using TransParmType = pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,
                                                          pele::physics::PhysicsType::transport_type>;

// what the advance of the lines of one level did, summed over the ranks
struct LineDiagnostics
{
  amrex::Long n_lines        = 0; // lines advanced
  amrex::Long n_skipped      = 0; // lines left untouched: unusable (no mass, rho or pressure)
  amrex::Long total_substeps = 0; // sub-steps over all the lines: the cost
  int max_substeps           = 0; // sub-steps of the worst line
  amrex::Long n_capped       = 0; // lines that hit setup::max_substeps: advanced ABOVE
                                  // their stability limit
  amrex::Long n_clipped      = 0; // elements where a Y_k < -clip_tol had to be floored
  amrex::Real max_neg_Y      = 0.0; // the worst negative Y_k floored
  // | sum m h after - before | / sum m |h| : at constant pressure the enthalpy
  // of a line, not its internal energy, is what diffusion must conserve
  amrex::Real drift_enthalpy = 0.0;
};

using ClemParticleContainer = amrex::AmrParticleContainer<RealData::ncomps, IntData::ncomps, 0, 0>;
using ClemParIter           = amrex::ParIter<RealData::ncomps, IntData::ncomps, 0, 0>;
using ClemParConstIter      = amrex::ParConstIter<RealData::ncomps, IntData::ncomps, 0, 0>;

class ClemContainer
{
public:
  using ParticleType = ClemParticleContainer::ParticleType;
  // per-grid fab of particle-index lists (indices into that grid's AoS)
  using CellVectors = std::map<int, amrex::BaseFab<std::vector<int>>>;

  ClemContainer() = default;

  // build the storage on the AMR hierarchy of amr
  void Define(amrex::AmrCore* amr);

  bool IsDefined() const { return m_pc != nullptr; }
  int FinestLevel() const { return m_pc->finestLevel(); }

  ClemParticleContainer& Particles() { return *m_pc; }
  const ClemParticleContainer& Particles() const { return *m_pc; }

  // ---- population of the lines ------------------------------------------------

  // InitializeFromLES: seed the lines of every cell of level lev from the LES
  // state `les` (on this level's particle BoxArray and DistributionMapping),
  // in conserved form: rho at rho_comp, T at temp_comp, rho e at eint_comp and
  // rho Y_k at spec_comp + k. Each cell gets setup::n_lem elements with
  //   rho = rho_LES,  vol = V_cell / n_lem,  mass = rho vol,
  //   T = T_LES,  eint = e_LES,  Y_k = Y_k,LES,  press = p(rho, T, Y),  phi = 0,
  //   lem_index = k,  entry_side = 0,
  // so the elements sum back to the cell. The level must be empty; the cell
  // map of lev is rebuilt on return.
  void InitializeFromLES(int lev, const amrex::MultiFab& les, int rho_comp, int temp_comp,
                         int eint_comp, int spec_comp);

  // ---- filter -------------------------------------------------------------------

  // GetFilteredScalars: write the filtered values of the lines of level lev
  // into `filtered` (RealData::ncomps components, indexed as RealData, on this
  // level's particle BoxArray and DistributionMapping). Per cell
  //   mass = sum_k m_k,   vol = sum_k v_k,   rho = mass / vol,
  //   q_tilde = sum_k m_k q_k / sum_k m_k   for q = T, eint, press, phi, Y_k.
  void GetFilteredScalars(int lev, amrex::MultiFab& filtered) const;

  // ---- cell map -----------------------------------------------------------------

  // (re)build the cell -> particle map of one level, or of all levels; must be
  // called after any Redistribute before the map is used
  void RebuildCellMap(int lev);
  void RebuildCellMaps();
  void InvalidateCellMaps();
  bool IsCellMapValid(int lev) const
  {
    return lev < static_cast<int>(m_cell_map_valid.size()) && m_cell_map_valid[lev] != 0;
  }

  // particle indices (into the AoS of grid grid_id on level lev) living in cell iv
  const std::vector<int>& ParticlesInCell(int lev, int grid_id, const amrex::IntVect& iv) const;

  // Redistribute of ONE level, Redistribute(lev, lev): the particles stay on
  // lev (delivered between its grids and ranks, the invalid ones removed). The
  // default AMReX Redistribute would move the lines of the covered cells to
  // the finer level. The map of lev is invalidated, so a stale map cannot be
  // used. Without lev: every level
  void Redistribute(int lev);
  void Redistribute();

  // ---- reactor of the lines (reactor/react.cpp) ------------------------------------

  // DiffusionStirringReactor: advance every line of level lev by the LES step
  // dt with the operators that act ALONG a line (the models of setup.H), in a
  // loop of explicit sub-steps sized by their stability limits. Nothing moves
  // between cells. trans_parm are the transport coefficients of PelePhysics.
  // The cell map of lev is rebuilt if stale.
  LineDiagnostics DiffusionStirringReactor(int lev, amrex::Real dt,
                                           const TransParmType* trans_parm);

  // ---- plotting (particle_plot.cpp) ----------------------------------------------

  // The particles go NEXT to the LES plotfile `grid_plotfile`, not inside it,
  // in a directory named after it (root + tag + step):
  //   LES plotfile        plt00048
  //   WritePlotFile       plt_particle_00048
  //   WriteLeafPlotFile   plt_particle_leave_00048
  // Inside, the layout ParaView's AMReX particle reader opens by default:
  //   Header               a copy of the header of the LES plotfile
  //   particles/Header     the particle header
  //   particles/Level_*/   the particle data
  // All the components are written, named as RealData / IntData, species as
  // Y(name). The LES plotfile must be complete (its Header written).

  // WritePlotFile: every particle of every level.
  void WritePlotFile(const std::string& grid_plotfile) const;

  // WriteLeafPlotFile: only the lines of the leaves: the particles of a cell
  // covered by the finer level are left out, so each point of the domain is
  // shown once, at its finest resolution.
  void WriteLeafPlotFile(const std::string& grid_plotfile);

private:
  // flip the sign of the id of every particle living in a covered cell: a
  // negative id hides the particle from the AMReX writers. Two calls restore
  // the container; nothing may redistribute in between
  void FlipCoveredParticles();

  std::unique_ptr<ClemParticleContainer> m_pc;
  // one map per AMR level, keyed by grid index
  amrex::Vector<CellVectors> m_cell_vectors;
  amrex::Vector<int> m_cell_map_valid; // per level; int, since Vector<bool> has no references
};

// slot k of the n of the line of cell iv: at fraction (k + 1) / (n + 1) of the
// cell along x, at its centre in the other directions; strictly interior, so
// AMReX bins it in its own cell. Sets the position and lem_index = k
inline void
PlaceOnLine(ClemContainer::ParticleType& p, const amrex::IntVect& iv, const int k, const int n,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cell_size)
{
  const amrex::Real dx_inner = cell_size[0] / static_cast<amrex::Real>(n + 1);
  AMREX_D_TERM(
    p.pos(0) = prob_lo[0] + iv[0] * cell_size[0] + static_cast<amrex::Real>(k + 1) * dx_inner;
    , p.pos(1) = prob_lo[1] + (iv[1] + 0.5) * cell_size[1];
    , p.pos(2) = prob_lo[2] + (iv[2] + 0.5) * cell_size[2];)
  p.idata(IntData::lem_index) = k;
}

} // namespace clem
#endif
