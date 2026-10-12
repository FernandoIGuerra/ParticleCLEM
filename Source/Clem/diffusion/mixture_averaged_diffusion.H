/**
 * This is a derived class from Diffusion.
 *
 * MixtureAveragedDiffusion: species and TEMPERATURE diffusion with the
 * mixture-averaged transport of PelePhysics, in mass coordinate and at
 * CONSTANT PRESSURE, on a PERIODIC line.
 *
 * With dm = rho A dx, d/dx = A rho d/dm:
 *
 *   DY_k/Dt = A^2 d/dm( S_k ) ,   S_k = F_k - Y_k sum_j F_j ,
 *   F_k     = ( rho Dhat_k ) dX_k/dm ,   Dhat_k = rho D_k W_k / W ,
 *
 *   DT/Dt   = ( A^2 / Cp ) d/dm( rho lambda dT/dm )      conduction
 *           + ( A^2 / Cp ) sum_k Cp_k S_k dT/dm          interspecies transport
 *
 * Dhat_k is what PelePhysics returns in Ddiag. sum_k S_k = 0, so the species
 * mass of the line is conserved and sum_k Y_k = 1 is kept. The conduction is
 * a flux divergence; the interspecies term is a gradient product, so the
 * temperature form does not conserve energy exactly: the enthalpy drift is
 * measured by the line advance (LineDiagnostics::drift_enthalpy).
 *
 * Isobaric closure: each sub-step holds the element pressure and, after T and
 * Y, recovers (rho, eint) from (p, T, Y). The mass of an element is fixed, so
 * its density, hence its volume, changes: the line breathes.
 *
 * Discretisation: central differences on the periodic faces, explicit Euler
 * in time. Face coefficients by harmonic mean, face Y_k and Cp_k by average.
 *
 * Stability: dt_sub = setup::diffusion_cfl dm_min^2 / max( beta_T, beta_Y ),
 *   beta_T = A^2 rho lambda / Cp ,   beta_Y = A^2 rho^2 D_k .
 */
#ifndef CLEM_MIXTURE_AVERAGED_DIFFUSION_H
#define CLEM_MIXTURE_AVERAGED_DIFFUSION_H

#include "diffusion.H"

namespace clem {

class MixtureAveragedDiffusion : public Diffusion<MixtureAveragedDiffusion>
{
  friend Diffusion;

protected:
  // X_k, Cp, Cp_k, W, rho Dhat_k, rho lambda of every element
  void UpdatePropertiesImplementation(LineScratch& line, const TransParmType* trans_parm) const;

  // the diffusion-number limit of the line; a huge value if nothing diffuses
  amrex::Real StableSubStepImplementation(const LineScratch& line, amrex::Real area2) const;

  // face fluxes S_k, conduction and interspecies product, then dY/dt, dT/dt.
  // Needs UpdateProperties first
  void ComputeRatesImplementation(LineScratch& line, amrex::Real area2) const;

  // T += dt dT/dt, Y += dt dY/dt (negative Y_k floored and renormalised),
  // then (rho, eint) from (p, T, Y)
  void AdvanceLineImplementation(LineScratch& line, amrex::Real dt_sub,
                                 LineDiagnostics& diag) const;

public:
  static constexpr const char* name = "mixture-averaged Y_k and T, mass coordinate, isobaric";

  MixtureAveragedDiffusion()                                           = default;
  ~MixtureAveragedDiffusion()                                          = default;
  MixtureAveragedDiffusion(MixtureAveragedDiffusion const&)            = delete;
  MixtureAveragedDiffusion& operator=(MixtureAveragedDiffusion const&) = delete;
  MixtureAveragedDiffusion(MixtureAveragedDiffusion&&)                 = delete;
  MixtureAveragedDiffusion& operator=(MixtureAveragedDiffusion&&)      = delete;
};

} // namespace clem
#endif
