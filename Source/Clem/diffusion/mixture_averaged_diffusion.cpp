#include <cmath>
#include <limits>

#include "mixture_averaged_diffusion.H"
#include "setup.H"

namespace clem {

namespace {
// face value of a transport coefficient: the series-resistance average, zero
// when either side vanishes
AMREX_FORCE_INLINE amrex::Real
HarmonicMean(const amrex::Real a, const amrex::Real b)
{
  const amrex::Real sum = a + b;
  return (sum > 0.0) ? (2.0 * a * b / sum) : 0.0;
}

// below this a negative Y_k is round-off: floored, but not counted as clipped
constexpr amrex::Real clip_tol = 1.0e-12;
} // namespace

void
MixtureAveragedDiffusion::UpdatePropertiesImplementation(LineScratch& line,
                                                         const TransParmType* trans_parm) const
{
  auto eos   = pele::physics::PhysicsType::eos();
  auto trans = pele::physics::PhysicsType::transport();

  // only the two coefficients the line equations need
  constexpr bool get_xi    = false;
  constexpr bool get_mu    = false;
  constexpr bool get_lam   = true;
  constexpr bool get_ddiag = true;
  constexpr bool get_chi   = false;

  for (int l = 0; l < line.n_elem; ++l) {
    const amrex::Real rho  = line.rho[l];
    const amrex::Real temp = line.T[l];
    const amrex::Real* Y   = &line.Y[l * NUM_SPECIES];

    eos.Y2X(Y, &line.X[l * NUM_SPECIES]);
    eos.Y2WBAR(Y, line.w_bar[l]);
    eos.RTY2Cp(rho, temp, Y, line.cp[l]);
    eos.T2Cpi(temp, &line.cp_k[l * NUM_SPECIES]);

    amrex::Real d_diag[NUM_SPECIES];
    amrex::Real mu       = 0.0;
    amrex::Real xi       = 0.0;
    amrex::Real lam      = 0.0;
    amrex::Real* chi_mix = nullptr;
    trans.transport(get_xi, get_mu, get_lam, get_ddiag, get_chi, temp, rho, Y, d_diag, chi_mix,
                    mu, xi, lam, trans_parm);

    // d/dx = A rho d/dm applied twice leaves one extra rho on each
    // coefficient; d_diag is already Dhat_k = rho D_k W_k / W
    line.rho_lambda[l] = rho * lam;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      line.rho_d[l * NUM_SPECIES + k] = rho * d_diag[k];
    }
  }
}

amrex::Real
MixtureAveragedDiffusion::StableSubStepImplementation(const LineScratch& line,
                                                      const amrex::Real area2) const
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real inv_mw[NUM_SPECIES];
  eos.inv_molecular_weight(inv_mw);

  // heat and species have their own diffusion number; the tighter one wins
  amrex::Real beta   = 0.0;
  amrex::Real dm_min = std::numeric_limits<amrex::Real>::max();

  for (int l = 0; l < line.n_elem; ++l) {
    dm_min = amrex::min<amrex::Real>(dm_min, line.mass[l]);

    // the T update goes through Cp: rho lambda / Cp
    if (line.cp[l] > 0.0) {
      beta = amrex::max<amrex::Real>(beta, area2 * line.rho_lambda[l] / line.cp[l]);
    }
    // the flux is driven by X_k but updates Y_k: rho Dhat_k W / W_k = rho^2 D_k
    for (int k = 0; k < NUM_SPECIES; ++k) {
      const amrex::Real rho2_d = line.rho_d[l * NUM_SPECIES + k] * line.w_bar[l] * inv_mw[k];
      beta                     = amrex::max<amrex::Real>(beta, area2 * rho2_d);
    }
  }

  if (beta <= 0.0) {
    return std::numeric_limits<amrex::Real>::max(); // no transport, no limit
  }
  return setup::diffusion_cfl * dm_min * dm_min / beta;
}

void
MixtureAveragedDiffusion::ComputeRatesImplementation(LineScratch& line,
                                                     const amrex::Real area2) const
{
  const int n_elem = line.n_elem;

  // periodic line: face f sits between the elements f and (f + 1) mod n_elem
  for (int f = 0; f < line.n_face; ++f) {
    const int lo = f;
    const int hi = (f + 1) % n_elem;
    // 1 / (mass-coordinate distance between the two element centres)
    const amrex::Real inv_dm = 2.0 / (line.mass[lo] + line.mass[hi]);

    amrex::Real* flux = &line.flux_species[f * NUM_SPECIES];

    // Fickian fluxes F_k = (rho Dhat_k) dX_k/dm and their sum, the correction
    // that makes the species fluxes close
    amrex::Real flux_sum = 0.0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      const amrex::Real coef  =  HarmonicMean(line.rho_d[lo * NUM_SPECIES + k], line.rho_d[hi * NUM_SPECIES + k]);
      flux[k]                 = coef * (line.X[hi * NUM_SPECIES + k] - line.X[lo * NUM_SPECIES + k]) * inv_dm;
      flux_sum                += flux[k];
    }

    // S_k = F_k - Y_k sum_j F_j, so sum_k S_k = 0; the interspecies transport
    // of the temperature equation, sum_k Cp_k S_k dT/dm, uses the same S_k
    const amrex::Real dT_dm   = (line.T[hi] - line.T[lo]) * inv_dm;
    amrex::Real interspecies = 0.0;
    
    for (int k = 0; k < NUM_SPECIES; ++k) {
      const amrex::Real Y_face  = 0.5 * (line.Y[lo * NUM_SPECIES + k] + line.Y[hi * NUM_SPECIES + k]);
      flux[k]                   -= Y_face * flux_sum;

      const amrex::Real cp_k_face   = 0.5 * (line.cp_k[lo * NUM_SPECIES + k] + line.cp_k[hi * NUM_SPECIES + k]);
      interspecies                  += cp_k_face * flux[k] * dT_dm;
    }

    line.flux_interspecies[f] = interspecies;

    // conduction only: the species enthalpy flux is the interspecies term
    line.flux_heat[f] = HarmonicMean(line.rho_lambda[lo], line.rho_lambda[hi]) * dT_dm;
  }

  // element l has the faces l (its l+1/2) and (l - 1) mod n_elem (its l-1/2)
  for (int l = 0; l < n_elem; ++l) {
    const int f_hi       = l;
    const int f_lo       = (l + n_elem - 1) % n_elem;
    const amrex::Real cp = (line.cp[l] > 0.0) ? line.cp[l] : 1.0;

    const amrex::Real conduction    = (line.flux_heat[f_hi] - line.flux_heat[f_lo]) / line.mass[l];
    const amrex::Real interspecies  = 0.5 * (line.flux_interspecies[f_hi] + line.flux_interspecies[f_lo]);
    
    line.dT_dt[l] = area2 * (conduction + interspecies) / cp;

    const amrex::Real fac = area2 / line.mass[l];
    
    for (int k = 0; k < NUM_SPECIES; ++k) {
      line.dY_dt[l * NUM_SPECIES + k] = fac * (line.flux_species[f_hi * NUM_SPECIES + k] - line.flux_species[f_lo * NUM_SPECIES + k]);
    }
  }
}

void
MixtureAveragedDiffusion::AdvanceLineImplementation(LineScratch& line, const amrex::Real dt_sub,
                                                    LineDiagnostics& diag) const
{
  auto eos = pele::physics::PhysicsType::eos();

  for (int l = 0; l < line.n_elem; ++l) {
    line.T[l] += dt_sub * line.dT_dt[l];

    amrex::Real* Y        = &line.Y[l * NUM_SPECIES];
    amrex::Real Y_sum     = 0.0;
    amrex::Real worst_neg = 0.0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      Y[k] += dt_sub * line.dY_dt[l * NUM_SPECIES + k];
      if (Y[k] < 0.0) {
        worst_neg = amrex::max<amrex::Real>(worst_neg, -Y[k]);
        Y[k]      = 0.0;
      }
      Y_sum += Y[k];
    }

    // sum_k Y_k = 1 holds analytically, so the renormalisation only undoes
    // what the flooring added
    if (worst_neg > 0.0) {
      const amrex::Real inv_sum = (Y_sum > 0.0) ? 1.0 / Y_sum : 0.0;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        Y[k] *= inv_sum;
      }
      diag.max_neg_Y = amrex::max<amrex::Real>(diag.max_neg_Y, worst_neg);
      if (worst_neg > clip_tol) {
        ++diag.n_clipped;
      }
    }

    // isobaric closure: rho and eint from the new (T, Y) at the held pressure
    eos.PYT2RE(line.press[l], Y, line.T[l], line.rho[l], line.eint[l]);
  }
}

} // namespace clem
