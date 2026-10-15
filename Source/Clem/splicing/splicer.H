/**
 * The splicing: the transport of the LEM lines between the cells of ONE level
 * by the mass flux of the LES step.
 *
 * Every cell c has faces with an outflux, N^+_c, and with an influx, N^-_c,
 * and each face its transfer mass of the step, mu_f = dt |mass_flux_f|. The
 * line of a cell is index-ordered: lem_index 0 is its head, the last its
 * tail. The splice is Sigma^omega = R_I^omega o A^omega:
 *
 *   A^omega    the donor serves its outflux faces in the order O^omega(N^+_c).
 *              Each face takes the segment of mass mu_f from the HEAD of what
 *              is left of the line; the element across the cut is split, both
 *              pieces keep their state, no average is formed. The segment goes
 *              to the neighbour through that face
 *   R_I^omega  the receptor keeps what is left of its line and appends the
 *              arriving blocks at its TAIL, ranked by O^omega(N^-_c); inside
 *              a block the order of the donor is kept
 *
 * O^omega is the policy of setup.H (reindexing/): Sigma^+, Sigma^- or the
 * two alternating.
 *
 * A face whose neighbour cell is not on the grids of the level (the boundary
 * of the domain, or a coarse/fine boundary) is a boundary of the lines:
 *
 *   outflux   the segment leaves the lines
 *   influx    one block of mass mu_f enters, with the state of the LES ghost
 *             cell behind the face (rho, T, Y_k; e and p from the EOS)
 *
 * The elements change in number and in mass; the regridding restores n_lem
 * equal-mass elements. Since the mass flux is the one of the LES density
 * update, after the splice the mass of every line is rho^{n+1} V_c.
 *
 * Functions:
 *
 * Splice          : the splice of the lines of a level
 *
 * LineMassDefect  : max | m_line - rho V_c | / rho V_c against a LES state
 *
 * LineTemperatureGap : max | mass-mean T of a line - T of its LES cell |
 */
#ifndef CLEM_SPLICER_H
#define CLEM_SPLICER_H

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>

#include "les_comps.H"
#include "particle_container.H"

namespace clem {

// what one Splice did, over the level and the ranks
struct SpliceDiagnostics
{
  amrex::Real mass_before = 0.0; // mass of the lines before
  amrex::Real mass_after  = 0.0; // and after
  amrex::Real mass_in     = 0.0; // entered through the boundary of the lines
  amrex::Real mass_out    = 0.0; // left through it
  amrex::Real max_export  = 0.0; // max over the cells of mu_out_c / m_c
  amrex::Long n_elements  = 0;   // elements on the level after
};

// Splice the lines of level lev over dt. mass_flux[d] is the mass flux of the
// step on the faces of direction d, signed along +d (PeleC::clem_mass_flux);
// les_border the LES state with at least one filled ghost cell; step the step
// count of the level, for the policy. Aborts if a cell exports more than it
// holds. On return the cell map of lev is valid and every line is relabelled
// lem_index = 0 .. N - 1, with N any number.
SpliceDiagnostics Splice(ClemContainer& container, int lev, amrex::Real dt,
                         const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& mass_flux,
                         const amrex::MultiFab& les_border, const LesComps& comps, int step);

// the largest relative difference between the mass of a line and the mass of
// its LES cell, rho V_c, over level lev
amrex::Real LineMassDefect(ClemContainer& container, int lev, const amrex::MultiFab& les,
                           int rho_comp);

// the largest difference between the mass-mean temperature of a line and the
// temperature of its LES cell, over level lev
amrex::Real LineTemperatureGap(ClemContainer& container, int lev, const amrex::MultiFab& les,
                               int temp_comp);

} // namespace clem
#endif
