/**
 * advance_clem.cpp: the step of the CLEM model inside PeleC.
 *
 *   PeleC::AdvanceClem  called by the MOL advance in place of
 *                       PeleC::react_state, with the same arguments
 *
 * The LES solves rho, rho u, rho E; the CLEM solves Y and T on the lines, in
 * a Lagrangian formulation. One step of a level is, in this order:
 *
 *      the halo of the level above is filled         (multiresolution/)
 *      from the lines of this level at the old time
 *   1. ForwardCoupling: isentropic pressure change   (coupling/)
 *   2. DiffusionStirringReactor: diffusion, reaction  (reactor/)
 *      (stirring to come) along the lines
 *   3. Splice: transport of the lines by the LES     (splicing/)
 *      mass flux, in the order of the policy
 *   4. Regrid: n_lem equal-mass elements per line    (regridding/)
 *   5. BackwardCoupling: conserved-delta write-back   (coupling/)
 *
 * The step ends with the regrid, so between steps every line of every level
 * holds n_lem equal-mass elements: what the line operators, the halo and the
 * average down assume.
 *
 * On entry (react_init = false) the MOL advance has left
 *   S_new   = U^n + dt ( F_AD + I_R )
 *   aux_src = F_AD, the advection-diffusion forcing of the step
 * and expects on return the new I_R (Reactions_Type) and
 *   S_new   = U^{n+1} = U^n + dt ( F_AD + I_R ).
 */
#include "PeleC.H"
#include "IndexDefines.H"

#include "kolmogorov_scale.H"
#include "stirring_rate.H"
#include "particle_container.H"
#include "setup.H"
#include "splicer.H"

void
PeleC::AdvanceClem(amrex::Real /*time*/, const amrex::Real dt, const bool react_init,
                   amrex::MultiFab* /*aux_src*/)
{
  BL_PROFILE("PeleC::AdvanceClem()");

  // No reaction on the LES: I_R = 0, so S_new = U^n + dt F_AD is already final
  if (react_init) {
    get_new_data(Reactions_Type).setVal(0.0);
    return;
  }

  amrex::MultiFab& S_new       = get_new_data(State_Type);
  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const clem::LesComps comps{URHO, UMX, UEDEN, UEINT, UTEMP, UFS};
  const clem::setup::CouplingType coupling;

  // The temperature of S_new is still that of the predictor
  computeTemp(S_new, 0);
  ProbeClem("LES after the MOL step");

  // The turbulent viscosity of the LES closure over this step, left by the
  // MOL stages in clem_nu_t (PeleC::getWALELESTerm)
  if (verbose != 0) {
    const auto n_cells = static_cast<amrex::Real>(clem_nu_t.boxArray().numPts());
    amrex::Print() << "CLEM: level " << level << " nu_T of the LES closure: max "
                   << clem_nu_t.max(0) << ", mean " << clem_nu_t.sum(0) / n_cells
                   << ", min " << clem_nu_t.min(0) << std::endl;
  }

  // The Kolmogorov scale of the cells, eta = (nu^3 / eps)^(1/4), from the
  // subgrid dissipation of the same closure (clem_eps) and the molecular
  // viscosity of S_new: the smallest eddy of the stirring
  const clem::KolmogorovDiagnostics kolmogorov =
    clem::KolmogorovScale(S_new, clem_eps, clem_eta, comps, geom.CellSize(0), trans_parms.device_parm());
  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " Kolmogorov scale: " << kolmogorov.n_turbulent
                   << " cells with eps > 0, eta min " << kolmogorov.eta_min << " max "
                   << kolmogorov.eta_max << ", Delta / eta max " << kolmogorov.delta_over_eta_max
                   << std::endl;
  }

  // The stirring of the lines: the event rate per unit length that gives the
  // D_T = nu_T / Sc_T of the closure, and the mean time between the events of
  // each line, dt_stir = 1 / (Lambda L)
  const clem::StirringDiagnostics stirring =
    clem::StirringRate(clem_nu_t, clem_eta, clem_lambda, clem_dt_stir, geom, dt);
  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " stirring: " << stirring.n_stirred
                   << " cells stirred (" << stirring.n_truncated
                   << " with eta below the line resolution), Lambda max " << stirring.lambda_max
                   << ", dt_stir min " << stirring.dt_stir_min << ", events per line in this step max "
                   << stirring.events_max << std::endl;
  }

  // The conserved filter of the lines before the CLEM step
  amrex::MultiFab conserved_pre;
  if (clem::setup::couple_back) {
    conserved_pre.define(grids, dmap, NUM_SPECIES + 1, 0);
    clem::ConservedSnapshot(*ClemPC, level, conserved_pre);
  }

  // The halo of the level above: a copy of the lines of this level at the old
  // time, before anything changes them
  if (level < parent->finestLevel()) {
    clem::HaloManager& halo = getLevel(level + 1).clem_halo;
    halo.FillFromCoarse(*ClemPC, level + 1);

    // and the level above starts the sum of its transfers over its substeps
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      getLevel(level + 1).clem_mass_transfer[dir].setVal(0.0);
    }

    if (verbose != 0) {
      amrex::Print() << "CLEM: level " << level + 1 << " halo filled from level " << level << ": "
                     << halo.NumHaloCells() << " cells, mass " << halo.TotalMass() << std::endl;
    }
  }

  // The mass through the faces of this level in this step, for the reflux of
  // the lines of the level below
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::MultiFab::Saxpy(clem_mass_transfer[dir], dt, clem_mass_flux[dir], 0, 0, 1, 0);
  }

  // 1. isentropic pressure change: the elements adopt the LES pressure
  coupling.ForwardCoupling(*ClemPC, level, S_new, comps);
  ProbeClem("after the forward coupling");

  // 2. diffusion, stirring and reaction along the lines of this level. The
  // change of [rho e]_c they make is what the LES energy equation does not
  // have (setup::energy_coupling = LesPlusSubgrid)
  amrex::MultiFab subgrid_energy(grids, dmap, 1, 0);
  amrex::MultiFab conserved_subgrid(grids, dmap, NUM_SPECIES + 1, 0);
  clem::ConservedSnapshot(*ClemPC, level, conserved_subgrid);
  amrex::MultiFab::Copy(subgrid_energy, conserved_subgrid, NUM_SPECIES, 0, 1, 0);

  const clem::LineDiagnostics diag =
    ClemPC->DiffusionStirringReactor(level, dt, trans_parms.device_parm(), clem_eta, clem_dt_stir);

  clem::ConservedSnapshot(*ClemPC, level, conserved_subgrid);
  amrex::MultiFab::LinComb(subgrid_energy, 1.0, conserved_subgrid, NUM_SPECIES, -1.0,
                           subgrid_energy, 0, 0, 1, 0);
  ProbeClem("after the line operators");

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " lines advanced for dt = " << dt << ": "
                   << diag.n_lines << " lines, sub-steps max " << diag.max_substeps << " mean "
                   << (diag.n_lines > 0 ? static_cast<amrex::Real>(diag.total_substeps) /
                                            static_cast<amrex::Real>(diag.n_lines)
                                        : 0.0)
                   << ", skipped " << diag.n_skipped << ", capped " << diag.n_capped
                   << ", clipped " << diag.n_clipped << " (max -Y " << diag.max_neg_Y
                   << "), enthalpy drift " << diag.drift_enthalpy << std::endl;
    if (clem::setup::react_lines) {
      amrex::Print() << "CLEM: level " << level << " reaction: " << diag.n_reacted
                     << " element sub-steps, RHS evaluations mean "
                     << (diag.n_reacted > 0 ? static_cast<amrex::Real>(diag.rhs_evals) /
                                                static_cast<amrex::Real>(diag.n_reacted)
                                            : 0.0)
                     << " max " << diag.max_rhs_evals << ", covered lines not reacted "
                     << diag.n_skip_react << std::endl;
    }
    if (clem::setup::stir_lines) {
      amrex::Print() << "CLEM: level " << level << " stirring: " << diag.stir_events
                     << " triplet maps on " << diag.n_stirred << " lines, max " << diag.max_stir_events
                     << " on one line" << std::endl;
    }
  }

  // 3. splicing: the lines follow the mass flux of this step. Sborder holds
  // the LES state with its ghost cells, the source of the inflow
  const clem::SpliceDiagnostics splice_diag = clem::Splice(
    *ClemPC, level, dt, clem_mass_flux, Sborder, comps, parent->levelSteps(level));
  ProbeClem("after the splice");

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " splice (" << clem::setup::ReindexerType::Name()
                   << "): mass " << splice_diag.mass_before << " -> " << splice_diag.mass_after
                   << ", in " << splice_diag.mass_in << ", out " << splice_diag.mass_out
                   << ", balance "
                   << (splice_diag.mass_after - splice_diag.mass_before - splice_diag.mass_in +
                       splice_diag.mass_out) /
                        splice_diag.mass_before
                   << ", max export " << splice_diag.max_export << ", elements "
                   << splice_diag.n_elements << std::endl;
  }

  // 4. regridding: every line back to n_lem elements of equal mass
  const clem::setup::RegridderType regridder;
  const clem::RegridDiagnostics regrid_diag = regridder.Regrid(*ClemPC, level);
  ProbeClem("after the regrid");

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " regrid: " << regrid_diag.n_added
                   << " elements added, " << regrid_diag.n_removed << " removed" << std::endl;
  }
  if (clem::setup::check_mass_flux) {
    amrex::Print() << "CLEM: level " << level << " line mass defect, max |m_line / (rho V) - 1| = "
                   << clem::LineMassDefect(*ClemPC, level, S_new, URHO) << std::endl;
  }

  // 5. conserved-delta write-back: the lines return rho Y_k and rho e
  amrex::MultiFab temp_before;
  if (verbose != 0 && clem::setup::couple_back) {
    temp_before.define(grids, dmap, 1, 0);
    amrex::MultiFab::Copy(temp_before, S_new, UTEMP, 0, 1, 0);
  }
  if (clem::setup::couple_back) {
    const clem::CouplingDiagnostics coupling_diag = coupling.BackwardCoupling(*ClemPC, level, S_new, S_old, conserved_pre, subgrid_energy, comps);

    if (verbose != 0) {
      amrex::Print() << "CLEM: level " << level << " write-back: max volume defect "
                     << coupling_diag.max_volume_defect << ", max projection "
                     << coupling_diag.max_projection << std::endl;
    }
  }
  computeTemp(S_new, 0);

  // 6. the lines take the energy of their cell
  if (clem::setup::couple_back && clem::setup::anchor_line_energy) {
    const clem::AnchorDiagnostics anchor = clem::AnchorLineEnergy(*ClemPC, level, S_new, comps);
    if (verbose != 0) {
      amrex::Print() << "CLEM: level " << level << " line energy anchored to the LES, max shift "
                     << anchor.max_shift << " K, energy moved " << anchor.energy_moved
                     << " erg (|.| " << anchor.energy_abs << ") of " << anchor.energy_line
                     << " erg, line mass " << anchor.mass << " g" << std::endl;
    }
  }
  ProbeClem("after the write-back");

  // the lines keep their own energy, so their temperature can drift from the
  // LES one: the largest gap between the mass-mean T of a line and its cell
  if (verbose != 0 && clem::setup::couple_back) {
    // how much the write-back moved the LES temperature, and where
    amrex::Real kick = 0.0;
    amrex::IntVect where(0);
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.validbox();
      const auto& before   = temp_before.const_array(mfi);
      const auto& after    = S_new.const_array(mfi, UTEMP);
      for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
        const amrex::Real d = std::abs(after(iv) - before(iv));
        if (!(d <= kick)) {
          kick  = std::isfinite(d) ? d : 1.0e300;
          where = iv;
        }
      }
    }
    amrex::Real kick_all = kick;
    amrex::ParallelDescriptor::ReduceRealMax(kick_all);
    if (kick == kick_all && kick > 0.0) {
      amrex::AllPrint() << "CLEM: level " << level << " write-back moved T_LES by up to " << kick
                        << " K at " << where << std::endl;
    }
    amrex::Print() << "CLEM: level " << level << " line vs LES temperature, max |T_line - T_LES| = "
                   << clem::LineTemperatureGap(*ClemPC, level, S_new, UTEMP) << " K" << std::endl;
  }
}
