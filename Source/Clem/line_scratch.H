/**
 * The workspace of the operators that act along one LEM line.
 *
 * A line operator cannot work on the particles in place: it needs the
 * elements of the cell in LINE order (the cell map gives them in storage
 * order), derived state that is not stored on the particle (X_k, Cp_k,
 * rho D_k, rho lambda, Cp) and face and rate buffers.
 *
 * LineScratch holds all of it. There is ONE per rank, grown on the first line
 * and reused for every line of every step: nothing is created or deleted
 * inside the loops.
 *
 * Layout: element quantities are indexed [l], species quantities
 * [l * NUM_SPECIES + k]. The line is PERIODIC: face f sits between the
 * elements f and (f + 1) mod n_elem, so there are n_elem faces.
 *
 * Functions:
 *
 * Get     : the workspace of this rank
 *
 * Gather  : copy the elements of one cell into the workspace, in line order
 *
 * Scatter : write the advanced state back to the particles
 */
#ifndef CLEM_LINE_SCRATCH_H
#define CLEM_LINE_SCRATCH_H

#include <algorithm>
#include <vector>

#include <AMReX_Vector.H>

#include "PelePhysics.H"

#include "particle_container.H"

namespace clem {

struct LineScratch
{
  using ParticleType = ClemContainer::ParticleType;

  // the workspace of this rank (of this thread, if the loop is ever threaded)
  static LineScratch& Get()
  {
    static thread_local LineScratch line;
    return line;
  }

  int n_elem = 0; // elements of the gathered line
  int n_face = 0; // faces = n_elem (periodic)

  // gathered element state, in line order (lem_index ascending)
  amrex::Vector<int> particle_index; // index of the element in the grid AoS
  amrex::Vector<amrex::Real> mass;   // element mass: the mass coordinate dm
  amrex::Vector<amrex::Real> rho;
  amrex::Vector<amrex::Real> T;
  amrex::Vector<amrex::Real> eint;
  amrex::Vector<amrex::Real> press; // held constant along the sub-steps
  amrex::Vector<amrex::Real> Y;

  // derived state, recomputed from (rho, T, Y) every sub-step
  amrex::Vector<amrex::Real> X;          // mole fractions
  amrex::Vector<amrex::Real> cp;         // mixture Cp
  amrex::Vector<amrex::Real> cp_k;       // species Cp_k
  amrex::Vector<amrex::Real> w_bar;      // mean molecular weight
  amrex::Vector<amrex::Real> rho_d;      // rho * (rho D_k W_k / W)
  amrex::Vector<amrex::Real> rho_lambda; // rho * lambda

  // face fluxes
  amrex::Vector<amrex::Real> flux_species;      // corrected species flux S_k
  amrex::Vector<amrex::Real> flux_heat;         // conduction (rho lambda) dT/dm
  amrex::Vector<amrex::Real> flux_interspecies; // sum_k Cp_k S_k dT/dm

  // rates of the current sub-step
  amrex::Vector<amrex::Real> dY_dt;
  amrex::Vector<amrex::Real> dT_dt;

  // make room for a line of n elements; after the first line this is only
  // the two assignments
  void Resize(const int n)
  {
    n_elem = n;
    n_face = n;
    if (static_cast<int>(mass.size()) >= n) {
      return;
    }
    const int n_spec = n * NUM_SPECIES;

    particle_index.resize(n);
    mass.resize(n);
    rho.resize(n);
    T.resize(n);
    eint.resize(n);
    press.resize(n);
    cp.resize(n);
    w_bar.resize(n);
    rho_lambda.resize(n);
    flux_heat.resize(n);
    flux_interspecies.resize(n);
    dT_dt.resize(n);

    Y.resize(n_spec);
    X.resize(n_spec);
    cp_k.resize(n_spec);
    rho_d.resize(n_spec);
    flux_species.resize(n_spec);
    dY_dt.resize(n_spec);
  }

  // copy the elements `in_cell` (indices into `particles`) into the
  // workspace, sorted by lem_index. Returns false, and the line must be
  // skipped, if it cannot be advanced: fewer than 2 elements, or an element
  // without mass, density or pressure
  bool Gather(const std::vector<int>& in_cell, const ParticleType* particles)
  {
    const int n = static_cast<int>(in_cell.size());
    if (n < 2) {
      return false;
    }
    Resize(n);

    std::copy(in_cell.begin(), in_cell.end(), particle_index.begin());
    std::sort(particle_index.begin(), particle_index.begin() + n, [&](const int a, const int b) {
      return particles[a].idata(IntData::lem_index) < particles[b].idata(IntData::lem_index);
    });

    for (int l = 0; l < n; ++l) {
      const auto& p = particles[particle_index[l]];
      if (p.rdata(RealData::mass) <= 0.0 || p.rdata(RealData::rho) <= 0.0 ||
          p.rdata(RealData::press) <= 0.0) {
        return false;
      }
      mass[l]  = p.rdata(RealData::mass);
      rho[l]   = p.rdata(RealData::rho);
      T[l]     = p.rdata(RealData::T);
      eint[l]  = p.rdata(RealData::eint);
      press[l] = p.rdata(RealData::press);
      for (int k = 0; k < NUM_SPECIES; ++k) {
        Y[l * NUM_SPECIES + k] = p.rdata(RealData::Y0 + k);
      }
    }
    return true;
  }

  // write (rho, vol, T, eint, Y) back. The mass and the pressure of an
  // element are not changed by the line operators; its density is, so its
  // volume follows as mass / rho
  void Scatter(ParticleType* particles) const
  {
    for (int l = 0; l < n_elem; ++l) {
      auto& p                = particles[particle_index[l]];
      p.rdata(RealData::rho)  = rho[l];
      p.rdata(RealData::vol)  = mass[l] / rho[l];
      p.rdata(RealData::T)    = T[l];
      p.rdata(RealData::eint) = eint[l];
      for (int k = 0; k < NUM_SPECIES; ++k) {
        p.rdata(RealData::Y0 + k) = Y[l * NUM_SPECIES + k];
      }
    }
  }

  // sum m h and sum m |h| of the gathered line, h = e + p / rho
  void Enthalpy(amrex::Real& sum, amrex::Real& sum_abs) const
  {
    for (int l = 0; l < n_elem; ++l) {
      const amrex::Real h = eint[l] + press[l] / rho[l];
      sum += mass[l] * h;
      sum_abs += mass[l] * std::abs(h);
    }
  }
};

} // namespace clem
#endif
