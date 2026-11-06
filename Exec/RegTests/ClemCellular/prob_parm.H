#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Vector.H>

// the columns of the ZND table after the distance behind the shock: rho, T,
// u_lab, then the NUM_SPECIES mass fractions
constexpr int ZND_RHO = 0;
constexpr int ZND_T = 1;
constexpr int ZND_U = 2;
constexpr int ZND_Y = 3;
constexpr int ZND_NCOL = 3 + NUM_SPECIES;

struct ProbParmDevice
{
  // the quiescent mixture ahead of the wave (CGS)
  amrex::Real p0 = 6.67e4; // 6.67 kPa
  amrex::Real T0 = 298.0;
  amrex::Real X_H2 = 2.0;
  amrex::Real X_O2 = 1.0;
  amrex::Real X_AR = 7.0;
  // the speed of the window, the frame of the simulation: the fresh mixture
  // enters at -D_window
  amrex::Real D_window = 1.61602e5; // D_CJ
  // the shock at the start
  amrex::Real x_shock = 4.5;
  // the random disturbance of the energy of the whole field at the start,
  // e* = e (1 + alpha f), f uniform on [-1, 1]
  amrex::Real alpha = 0.1;
  // an optional hot spot behind the shock at the start: a disc of radius
  // hs_r, its centre hs_x behind the shock and at height hs_y, raised to the
  // pressure hs_p at the density of the ZND state (hs_p <= 0: none)
  amrex::Real hs_p = 0.0;
  amrex::Real hs_x = 0.1;
  amrex::Real hs_y = 0.3;
  amrex::Real hs_r = 0.1;
  // filled in amrex_probinit
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac0 = {0.0};
  amrex::Real rho0 = 0.0;
  amrex::Real e0 = 0.0;
  int znd_n = 0;
  amrex::Real* d_znd_x = nullptr;
  amrex::Real* d_znd_q = nullptr; // [ZND_NCOL * i + column]
};

struct ProbParmHost
{
  std::string profile = "znd_davis_2H2_O2_7Ar.dat";
  amrex::Vector<amrex::Real> h_znd_x;
  amrex::Vector<amrex::Real> h_znd_q;
  amrex::Gpu::DeviceVector<amrex::Real> znd_x;
  amrex::Gpu::DeviceVector<amrex::Real> znd_q;
  ProbParmHost() : znd_x(0), znd_q(0) {}
};
#endif
