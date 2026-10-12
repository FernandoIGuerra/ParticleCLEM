#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

// Hydrogen jet in supersonic crossflow (JISCF) in a rectangular combustor.
// Reproduces the configuration of Jin, Cai, Hong, Zhang & Liang,
// Combustion and Flame 268 (2024) 113627.
//
// All quantities are CGS, as everywhere else in PeleC:
//   length  cm      time  s        mass  g
//   density g/cm^3  pres  dyn/cm^2 (1 atm = 1.01325e6)  temperature K
struct ProbParmDevice
{
  // ---- Vitiated-air crossflow entering at xlo (Table 1, Case 1) ----
  amrex::Real p_in = 1.01325e6; // 1 atm
  amrex::Real T_in = 1000.0;    // K
  amrex::Real u_in = 9.30e4;    // 930 m/s -> M = 1.5
  amrex::Real YO2_in = 0.23;
  amrex::Real YN2_in = 0.77;

  // ---- Transverse H2 jet issuing from the bottom (ylo) wall ----
  amrex::Real p_jet = 1.30709e6; // 1.29 atm
  amrex::Real T_jet = 300.0;     // K
  amrex::Real v_jet = 1.321e5;   // 1321 m/s -> M = 1.02
  amrex::Real jet_x = 4.2;       // injector centre, 42 mm from the inlet
  amrex::Real jet_z = 0.25;      // injector centre, spanwise mid-plane
  amrex::Real jet_d = 0.1;       // injector diameter, 1 mm

  // Injector start-up. jet_ramp > 0 ramps the injection velocity linearly
  // from 0 to v_jet over that interval, which softens the start-up
  // transient. Defaults reproduce the paper: the jet is on from t = 0.
  amrex::Real jet_start = 0.0;
  amrex::Real jet_ramp = 0.0;

  // ---- Derived in amrex_probinit, not read from the inputs file ----
  amrex::Real rho_in = 0.0;
  amrex::Real e_in = 0.0;
  amrex::Real rho_jet = 0.0;
  amrex::Real e_jet = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_in = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_jet = {{0.0}};
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
