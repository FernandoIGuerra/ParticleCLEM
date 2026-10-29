#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  // the mixture, as mole fractions (normalised in amrex_probinit)
  amrex::Real X_H2 = 2.0;
  amrex::Real X_O2 = 1.0;
  amrex::Real X_AR = 7.0;
  // the quiescent mixture ahead of the detonation (CGS)
  amrex::Real p0 = 6.67e4; // 6.67 kPa
  amrex::Real T0 = 298.0;
  // the driver at the closed end, x < x_driver: the same mixture, hot and at
  // high pressure, which ignites and launches the detonation
  amrex::Real p_driver = 1.0e6;
  amrex::Real T_driver = 2500.0;
  amrex::Real x_driver = 0.5;
  // filled in amrex_probinit
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  amrex::Real rho0 = 0.0;
  amrex::Real e0 = 0.0;
  amrex::Real rho_driver = 0.0;
  amrex::Real e_driver = 0.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
