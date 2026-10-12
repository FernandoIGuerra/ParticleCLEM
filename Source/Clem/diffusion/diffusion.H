/**
 * Here we create the base class of the molecular diffusion along a LEM line.
 *
 * Diffusion<Derived>: the base of every diffusion model (CRTP, no virtual
 * calls). A model moves species mass and heat BETWEEN the elements of one
 * line, never between cells (that is the splicing). It works on the gathered
 * line (LineScratch) and is driven by the sub-step loop of the line advance
 * (advance_clem.cpp), one explicit sub-step at a time:
 *
 *   UpdateProperties  the coefficients and driving forces from (rho, T, Y)
 *   StableSubStep     the largest stable sub-step of the line
 *   ComputeRates      the face fluxes, then dY/dt and dT/dt
 *   AdvanceLine       T, Y advanced by dt_sub; rho and eint recovered
 *
 * area2 is the square of the cross-section of the line, A = V_cell / length.
 *
 * The model is chosen at compile time in setup.H.
 */
#ifndef CLEM_DIFFUSION_H
#define CLEM_DIFFUSION_H

#include <AMReX_REAL.H>

#include "line_scratch.H"

namespace clem {

template <typename DerivedDiffusion>
class Diffusion
{
  friend DerivedDiffusion;

  // only a derived model can build the base
  Diffusion() = default;

  const DerivedDiffusion& derived() const { return static_cast<const DerivedDiffusion&>(*this); }

public:
  ~Diffusion()                           = default;
  Diffusion(Diffusion const&)            = delete;
  Diffusion& operator=(Diffusion const&) = delete;
  Diffusion(Diffusion&&)                 = delete;
  Diffusion& operator=(Diffusion&&)      = delete;

  void UpdateProperties(LineScratch& line, const TransParmType* trans_parm) const
  {
    derived().UpdatePropertiesImplementation(line, trans_parm);
  }

  amrex::Real StableSubStep(const LineScratch& line, const amrex::Real area2) const
  {
    return derived().StableSubStepImplementation(line, area2);
  }

  void ComputeRates(LineScratch& line, const amrex::Real area2) const
  {
    derived().ComputeRatesImplementation(line, area2);
  }

  void AdvanceLine(LineScratch& line, const amrex::Real dt_sub, LineDiagnostics& diag) const
  {
    derived().AdvanceLineImplementation(line, dt_sub, diag);
  }

  static constexpr const char* Name() { return DerivedDiffusion::name; }
};

} // namespace clem
#endif
