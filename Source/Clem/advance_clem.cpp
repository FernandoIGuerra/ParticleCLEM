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
 *   2. DiffusionStirringReactor                       (reactor/react.cpp)
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

    if (verbose != 0) {
      amrex::Print() << "CLEM: level " << level + 1 << " halo filled from level " << level << ": "
                     << halo.NumHaloCells() << " cells, mass " << halo.TotalMass() << std::endl;
    }
  }

  // 1. isentropic pressure change: the elements adopt the LES pressure
  coupling.ForwardCoupling(*ClemPC, level, S_new, comps);

  // 2. diffusion, stirring and reaction along the lines of this level
  const clem::LineDiagnostics diag =
    ClemPC->DiffusionStirringReactor(level, dt, trans_parms.device_parm());

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " lines advanced for dt = " << dt << ": "
                   << diag.n_lines << " lines, sub-steps max " << diag.max_substeps << " mean "
                   << (diag.n_lines > 0 ? static_cast<amrex::Real>(diag.total_substeps) /
                                            static_cast<amrex::Real>(diag.n_lines)
                                        : 0.0)
                   << ", skipped " << diag.n_skipped << ", capped " << diag.n_capped
                   << ", clipped " << diag.n_clipped << " (max -Y " << diag.max_neg_Y
                   << "), enthalpy drift " << diag.drift_enthalpy << std::endl;
  }

  // 3. splicing: the lines follow the mass flux of this step. Sborder holds
  // the LES state with its ghost cells, the source of the inflow
  const clem::SpliceDiagnostics splice_diag = clem::Splice(
    *ClemPC, level, dt, clem_mass_flux, Sborder, comps, parent->levelSteps(level));

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

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " regrid: " << regrid_diag.n_added
                   << " elements added, " << regrid_diag.n_removed << " removed" << std::endl;
  }
  if (clem::setup::check_mass_flux) {
    amrex::Print() << "CLEM: level " << level << " line mass defect, max |m_line / (rho V) - 1| = "
                   << clem::LineMassDefect(*ClemPC, level, S_new, URHO) << std::endl;
  }

  // 5. conserved-delta write-back: the lines return rho Y_k and rho e
  if (clem::setup::couple_back) {
    const clem::CouplingDiagnostics coupling_diag = coupling.BackwardCoupling(*ClemPC, level, S_new, S_old, conserved_pre, comps);

    if (verbose != 0) {
      amrex::Print() << "CLEM: level " << level << " write-back: max volume defect "
                     << coupling_diag.max_volume_defect << ", max projection "
                     << coupling_diag.max_projection << std::endl;
    }
  }
}
