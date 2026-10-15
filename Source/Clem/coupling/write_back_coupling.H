/**
 * This is a derived class from Coupling.
 *
 * WriteBackCoupling: isentropic pressure forward, conserved-delta write-back.
 *
 * ForwardCoupling: every element adopts the pressure P of its LES cell, from
 * ITS OWN pressure p, as isentropic work at frozen composition, with the
 * local gamma = Cp / Cv of the element:
 *
 *   T   <- T   ( P / p )^( (gamma - 1) / gamma )
 *   rho <- rho ( P / p )^( 1 / gamma ) ,   vol = mass / rho ,   p <- P
 *   e   += e( rho_new, T_new, Y ) - e( rho_old, T_old, Y )
 *
 * The energy is updated in difference form, so a constant bias of the
 * temperature representation cannot accumulate.
 *
 * BackwardCoupling: the LES keeps rho and rho u; the lines return rho Y_k
 * and rho e as the change of their conserved filter over the CLEM step:
 *
 *   1. volume renormalisation: the lines breathe at constant pressure, so
 *      sum_p v_p drifts from V_c. Every element of the cell is rescaled by
 *      the same isentropic ratio f = V_c / sum_p v_p,
 *        rho <- rho / f ,  T <- T f^-(gamma - 1) ,  p <- p f^-gamma ,
 *        e   += e( new ) - e( old ) ,
 *      which turns the p dV work of the breathing back into internal energy
 *   2. Delta[rho phi]_c = [rho phi]_c^post - [rho phi]_c^pre ,  phi = Y_k, e
 *   3. ( rho Y_k )^{n+1} = ( rho Y_k )^n + Delta[rho Y_k]_c , floored at 0 and
 *      PROJECTED onto the LES density, sum_k ( rho Y_k )^{n+1} = rho^{n+1}
 *   4. ( rho e )^{n+1} = ( rho e )^n + Delta[rho e]_c , and rho E rebuilt
 *      with the kinetic energy of the LES
 *
 * The delta carries everything the lines did in the step: isentropic work,
 * diffusion, stirring, reaction and the spliced transport.
 */
#ifndef CLEM_WRITE_BACK_COUPLING_H
#define CLEM_WRITE_BACK_COUPLING_H

#include "coupling.H"

namespace clem {

class WriteBackCoupling : public Coupling<WriteBackCoupling>
{
  friend Coupling;

protected:
  void ForwardCouplingImplementation(ClemContainer& container, int lev,
                                     const amrex::MultiFab& les_new, const LesComps& comps) const;

  CouplingDiagnostics BackwardCouplingImplementation(ClemContainer& container, int lev,
                                                     amrex::MultiFab& les_new,
                                                     const amrex::MultiFab& les_old,
                                                     const amrex::MultiFab& conserved_pre,
                                                     const amrex::MultiFab& subgrid_energy,
                                                     const LesComps& comps) const;

private:
  // step 1 of BackwardCoupling; returns max | f - 1 | over the level (local
  // to the rank)
  amrex::Real RenormalizeVolume(ClemContainer& container, int lev) const;

public:
  static constexpr const char* name = "isentropic pressure forward, conserved-delta write-back";

  WriteBackCoupling()                                    = default;
  ~WriteBackCoupling()                                   = default;
  WriteBackCoupling(WriteBackCoupling const&)            = delete;
  WriteBackCoupling& operator=(WriteBackCoupling const&) = delete;
  WriteBackCoupling(WriteBackCoupling&&)                 = delete;
  WriteBackCoupling& operator=(WriteBackCoupling&&)      = delete;
};

} // namespace clem
#endif
