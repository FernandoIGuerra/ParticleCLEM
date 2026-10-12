/**
 * Here we create the base class of the coupling between the LES and the lines.
 *
 * Coupling<Derived>: the base of every coupling strategy (CRTP, no virtual
 * calls). The LES solves rho, rho u, rho E; the lines solve Y and T. A
 * strategy owns the two directions, on the lines of ONE level:
 *
 *   ForwardCoupling   LES -> lines, at the start of the CLEM step: what the
 *                     resolved field imposes on every element
 *   BackwardCoupling  lines -> LES, at the end of the CLEM step: what the
 *                     lines return to the LES state
 *
 * and both use the conserved filter of a line,
 *
 *   [rho phi]_c = sum_p m_p phi_p / V_c ,
 *
 * given by ConservedSnapshot for phi = Y_k and e.
 *
 * The strategy is chosen at compile time in setup.H.
 */
#ifndef CLEM_COUPLING_H
#define CLEM_COUPLING_H

#include <AMReX_MultiFab.H>

#include "les_comps.H"
#include "particle_container.H"

namespace clem {

// what one BackwardCoupling did, over the level and the ranks
struct CouplingDiagnostics
{
  amrex::Real max_volume_defect = 0.0; // max | V_c / sum_p v_p - 1 | found
  amrex::Real max_projection    = 0.0; // max | rho_LES / sum_k (rho Y_k) - 1 | applied
};

// the conserved filter of the lines of level lev: conserved(.., k) =
// [rho Y_k]_c for k < NUM_SPECIES and conserved(.., NUM_SPECIES) = [rho e]_c.
// conserved has NUM_SPECIES + 1 components on the particle grids of lev
void ConservedSnapshot(const ClemContainer& container, int lev, amrex::MultiFab& conserved);

template <typename DerivedCoupling>
class Coupling
{
  friend DerivedCoupling;

  // only a derived strategy can build the base
  Coupling() = default;

  const DerivedCoupling& derived() const { return static_cast<const DerivedCoupling&>(*this); }

public:
  ~Coupling()                          = default;
  Coupling(Coupling const&)            = delete;
  Coupling& operator=(Coupling const&) = delete;
  Coupling(Coupling&&)                 = delete;
  Coupling& operator=(Coupling&&)      = delete;

  // LES -> lines. les_new is the LES state at the end of the LES step, with a
  // temperature consistent with it
  void ForwardCoupling(ClemContainer& container, const int lev, const amrex::MultiFab& les_new,
                       const LesComps& comps) const
  {
    derived().ForwardCouplingImplementation(container, lev, les_new, comps);
  }

  // lines -> LES. conserved_pre is the ConservedSnapshot taken at the start
  // of the CLEM step, before ForwardCoupling; les_old the LES state at the
  // start of the LES step. les_new is modified
  CouplingDiagnostics BackwardCoupling(ClemContainer& container, const int lev,
                                       amrex::MultiFab& les_new, const amrex::MultiFab& les_old,
                                       const amrex::MultiFab& conserved_pre,
                                       const LesComps& comps) const
  {
    return derived().BackwardCouplingImplementation(container, lev, les_new, les_old,
                                                    conserved_pre, comps);
  }

  static constexpr const char* Name() { return DerivedCoupling::name; }
};

} // namespace clem
#endif
