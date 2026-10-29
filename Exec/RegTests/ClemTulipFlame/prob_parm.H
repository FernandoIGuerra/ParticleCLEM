#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  // H2/air of equivalence ratio phi, at rest (CGS)
  amrex::Real phi = 1.0;
  amrex::Real p0 = 1.01325e6; // 1 atm
  amrex::Real T0 = 300.0;
  // the ignition kernel: the same mixture at T_ign in a cylinder (along z) of
  // radius r_ign centred on the closed end x = 0, mid-height
  amrex::Real T_ign = 1800.0;
  amrex::Real r_ign = 0.05;
  // filled in amrex_probinit
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  amrex::Real rho0 = 0.0;
  amrex::Real e0 = 0.0;
  amrex::Real rho_ign = 0.0;
  amrex::Real e_ign = 0.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
