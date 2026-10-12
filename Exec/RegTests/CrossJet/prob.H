#ifndef PROB_H
#define PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

// Uniform crossflow everywhere. The paper states that the initial pressure,
// temperature, Mach number and species coincide with the inflow conditions,
// so the wall boundary layers grow from a slug profile as the run proceeds.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real rho = prob_parm.rho_in;
  const amrex::Real u = prob_parm.u_in;

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * prob_parm.e_in;
  state(i, j, k, UEDEN) = rho * (prob_parm.e_in + 0.5 * u * u);
  state(i, j, k, UTEMP) = prob_parm.T_in;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * prob_parm.Y_in[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* x,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  // On entry s_ext holds a reflection of the interior state (odd for the
  // momenta, even for everything else), i.e. an adiabatic no-slip wall.
  // Anything left untouched below stays a wall.

  if ((idir == 0) && (sgn == 1)) {
    // xlo: vitiated-air crossflow. The inflow is supersonic, so every
    // characteristic points into the domain and the full state is imposed.
    const amrex::Real rho = prob_parm.rho_in;
    const amrex::Real u = prob_parm.u_in;

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = 0.0;
    s_ext[UMZ] = 0.0;
    s_ext[UEINT] = rho * prob_parm.e_in;
    s_ext[UEDEN] = rho * (prob_parm.e_in + 0.5 * u * u);
    s_ext[UTEMP] = prob_parm.T_in;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] = rho * prob_parm.Y_in[n];
    }

  } else if ((idir == 1) && (sgn == 1)) {
    // ylo: adiabatic no-slip wall pierced by the circular H2 injector.
    // In 2D the injector degenerates to a slot of width jet_d.
    const amrex::Real dxj = x[0] - prob_parm.jet_x;
#if AMREX_SPACEDIM == 3
    const amrex::Real dzj = x[2] - prob_parm.jet_z;
    const amrex::Real r2 = dxj * dxj + dzj * dzj;
#else
    const amrex::Real r2 = dxj * dxj;
#endif
    const amrex::Real rj = 0.5 * prob_parm.jet_d;

    if ((r2 < rj * rj) && (time >= prob_parm.jet_start)) {
      amrex::Real ramp = 1.0;
      if (prob_parm.jet_ramp > 0.0) {
        ramp = amrex::min<amrex::Real>( 1.0, (time - prob_parm.jet_start) / prob_parm.jet_ramp);
      }
      // The jet is sonic/supersonic, so the full state is imposed here too.
      const amrex::Real rho = prob_parm.rho_jet;
      const amrex::Real v = ramp * prob_parm.v_jet;

      s_ext[URHO] = rho;
      s_ext[UMX] = 0.0;
      s_ext[UMY] = rho * v;
      s_ext[UMZ] = 0.0;
      s_ext[UEINT] = rho * prob_parm.e_jet;
      s_ext[UEDEN] = rho * (prob_parm.e_jet + 0.5 * v * v);
      s_ext[UTEMP] = prob_parm.T_jet;
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[UFS + n] = rho * prob_parm.Y_jet[n];
      }
    }
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
