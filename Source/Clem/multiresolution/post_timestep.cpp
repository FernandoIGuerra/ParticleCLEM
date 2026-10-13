/**
 * post_timestep.cpp: the coarse/fine sync of the LEM lines inside PeleC.
 *
 *   PeleC::PostTimestepClem  called by PeleC::post_timestep of a level with a
 *                            finer level above it, once both are at the same
 *                            time and the LES is refluxed and averaged down
 *
 * In this order:
 *
 *   1. RefluxMass: the lines of the cells next to     (refluxer.H)
 *      the fine level get the mass the fine level
 *      exchanged, not the one their own splice did
 *   2. Regrid of those cells only                     (regridding/)
 *   3. conserved-delta write-back in those cells:
 *        rho Y_k += Delta[rho Y_k]_reflux, projected onto the LES density
 *        rho e   += Delta[rho e]_reflux
 *      The LES reflux no longer carries the species (Diffusion.cpp), so this
 *      is their only correction. rho E is the one PeleC refluxed.
 *   4. AverageDown: the covered cells of this level   (averager.H)
 *      from the lines of the level above
 */
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include "PeleC.H"
#include "IndexDefines.H"

#include "averager.H"
#include "coupling.H"
#include "refluxer.H"
#include "setup.H"
#include "splicer.H"

void
PeleC::PostTimestepClem()
{
  BL_PROFILE("PeleC::PostTimestepClem()");
  if (ClemPC == nullptr || !ClemPC->IsDefined() || level >= parent->finestLevel()) {
    return;
  }

  PeleC& fine                    = getLevel(level + 1);
  const amrex::IntVect ref_ratio = parent->refRatio(level);
  amrex::MultiFab& S_new         = get_new_data(State_Type);

  // The transfer still to be made through the interface faces: what the fine
  // level moved over its substeps minus what this level moved in its step
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> mismatch;
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    mismatch[dir].define(getEdgeBoxArray(dir), dmap, 1, 0);
    mismatch[dir].setVal(0.0);
  }
  amrex::average_down_faces(
    amrex::GetArrOfConstPtrs(fine.clem_mass_transfer), amrex::GetArrOfPtrs(mismatch), ref_ratio,
    geom);
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    // the mean over the sub-faces of a face, times their number, is their sum
    amrex::Real n_sub_faces = 1.0;
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
      n_sub_faces *= (d == dir) ? 1.0 : static_cast<amrex::Real>(ref_ratio[d]);
    }
    mismatch[dir].mult(n_sub_faces);
    amrex::MultiFab::Saxpy(mismatch[dir], -parent->dtLevel(level), clem_mass_flux[dir], 0, 0, 1, 0);
  }

  // 1 in the cells covered by the level above, with a ghost cell
  const amrex::iMultiFab covered = amrex::makeFineMask(
    grids, dmap, amrex::IntVect(1), fine.boxArray(), ref_ratio, geom.periodicity(), 0, 1);

  amrex::MultiFab conserved_pre(grids, dmap, NUM_SPECIES + 1, 0);
  clem::ConservedSnapshot(*ClemPC, level, conserved_pre);

  // 1. the mass reflux of the lines
  amrex::iMultiFab touched(grids, dmap, 1, 0);
  const clem::RefluxDiagnostics reflux_diag =
    clem::RefluxMass(*ClemPC, level, mismatch, covered, parent->levelSteps(level), touched);

  // 2. the lines that changed back to n_lem equal-mass elements
  const clem::setup::RegridderType regridder;
  regridder.Regrid(*ClemPC, level, &touched);

  // 3. what the reflux did to the conserved filter of those lines goes to the LES
  amrex::MultiFab conserved_post(grids, dmap, NUM_SPECIES + 1, 0);
  clem::ConservedSnapshot(*ClemPC, level, conserved_post);

  amrex::Real max_projection = 0.0;
  for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const auto& state    = S_new.array(mfi);
    const auto& pre      = conserved_pre.const_array(mfi);
    const auto& post     = conserved_post.const_array(mfi);
    const auto& mark     = touched.const_array(mfi);

    amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
      if (mark(i, j, k) != 1) {
        return;
      }
      amrex::Real rho_Y_sum = 0.0;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, UFS + n) = amrex::max<amrex::Real>(
          0.0, state(i, j, k, UFS + n) + (post(i, j, k, n) - pre(i, j, k, n)));
        rho_Y_sum += state(i, j, k, UFS + n);
      }
      const amrex::Real projection = state(i, j, k, URHO) / rho_Y_sum;
      max_projection = amrex::max<amrex::Real>(max_projection, std::abs(projection - 1.0));
      for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i, j, k, UFS + n) *= projection;
      }
      state(i, j, k, UEINT) += post(i, j, k, NUM_SPECIES) - pre(i, j, k, NUM_SPECIES);
    });
  }
  amrex::ParallelDescriptor::ReduceRealMax(max_projection);

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " mass reflux: " << reflux_diag.n_cells
                   << " cells, removed " << reflux_diag.mass_removed << ", added "
                   << reflux_diag.mass_added << ", max |dm| / m " << reflux_diag.max_fraction
                   << ", max projection " << max_projection << std::endl;
  }

  // 4. the covered cells take the average of their children
  const clem::AverageDiagnostics diag = clem::AverageDown(*ClemPC, level);

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " average down from level " << level + 1 << ": "
                   << diag.n_covered << " covered cells, " << diag.n_rebuilt << " rebuilt"
                   << std::endl;
  }
  if (clem::setup::check_mass_flux) {
    amrex::Print() << "CLEM: level " << level
                   << " line mass defect after the sync, max |m_line / (rho V) - 1| = "
                   << clem::LineMassDefect(*ClemPC, level, S_new, URHO) << std::endl;
  }
}
