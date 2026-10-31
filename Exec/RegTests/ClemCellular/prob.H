#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
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

// the fresh mixture of the frame of the window, moving at -D_window
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
FreshState(amrex::Real* s, ProbParmDevice const& prob_parm)
{
  const amrex::Real rho = prob_parm.rho0;
  const amrex::Real u = -prob_parm.D_window;
  s[URHO] = rho;
  s[UMX] = rho * u;
  s[UMY] = 0.0;
  s[UMZ] = 0.0;
  s[UEINT] = rho * prob_parm.e0;
  s[UEDEN] = rho * (prob_parm.e0 + 0.5 * u * u);
  s[UTEMP] = prob_parm.T0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s[UFS + n] = rho * prob_parm.massfrac0[n];
  }
}

// a uniform number in [-1, 1] from the cell index: the same on any number of
// ranks, and different from cell to cell
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
CellNoise(const int i, const int j)
{
  auto h = static_cast<unsigned int>(i) * 73856093U ^
           static_cast<unsigned int>(j) * 19349663U;
  h ^= h >> 13;
  h *= 0x5bd1e995U;
  h ^= h >> 15;
  return 2.0 * static_cast<amrex::Real>(h) / 4294967295.0 - 1.0;
}

// the 1-D ZND detonation in the frame of the window, the shock at x_shock and
// the fresh mixture ahead, with random disturbances of the energy behind the
// shock, which the transverse waves grow from
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real behind = prob_parm.x_shock - x;

  if (behind <= 0.0) {
    amrex::Real s[NVAR] = {0.0};
    FreshState(s, prob_parm);
    for (int n = 0; n < NVAR; n++) {
      state(i, j, k, n) = s[n];
    }
    return;
  }

  // the table is ascending in the distance behind the shock
  const int n_table = prob_parm.znd_n;
  const amrex::Real* xt = prob_parm.d_znd_x;
  const amrex::Real* q = prob_parm.d_znd_q;
  int lo = 0;
  int hi = n_table - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (xt[mid] <= behind) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  const amrex::Real w =
    amrex::max<amrex::Real>(amrex::min<amrex::Real>((behind - xt[lo]) / (xt[hi] - xt[lo]), 1.0), 0.0);
  const auto at = [=](const int c) {
    return (1.0 - w) * q[ZND_NCOL * lo + c] + w * q[ZND_NCOL * hi + c];
  };

  const amrex::Real rho = at(ZND_RHO);
  const amrex::Real u = at(ZND_U) - prob_parm.D_window;
  amrex::Real T = at(ZND_T);
  amrex::Real massfrac[NUM_SPECIES];
  amrex::Real sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = amrex::max<amrex::Real>(at(ZND_Y + n), 0.0);
    sum += massfrac[n];
  }
  for (amrex::Real& Yn : massfrac) {
    Yn /= sum;
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real e = 0.0;
  eos.TY2E(T, massfrac, e);
  e *= 1.0 + prob_parm.alpha * CellNoise(i, j);
  eos.REY2T(rho, e, massfrac, T);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * e;
  state(i, j, k, UEDEN) = rho * (e + 0.5 * u * u);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

// the fresh mixture enters supersonically through the high x side
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* /*s_int[NVAR]*/,
  amrex::Real* s_ext,
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  if ((idir == 0) && (sgn == -1)) {
    FreshState(s_ext, prob_parm);
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
