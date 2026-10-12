/**
 * post_timestep.cpp: the coarse/fine sync of the LEM lines inside PeleC.
 *
 *   PeleC::PostTimestepClem  called by PeleC::post_timestep of a level with a
 *                            finer level above it, once both are at the same
 *                            time and the LES is refluxed and averaged down
 *
 * In this order:
 *
 *   1. reflux of the lines next to the fine level     [to come]
 *   2. AverageDown: the covered cells of this level   (averager.H)
 *      from the lines of the level above
 */
#include "PeleC.H"

#include "averager.H"

void
PeleC::PostTimestepClem()
{
  BL_PROFILE("PeleC::PostTimestepClem()");
  if (ClemPC == nullptr || !ClemPC->IsDefined() || level >= parent->finestLevel()) {
    return;
  }

  // 1. reflux of the lines

  // 2. the covered cells take the average of their children
  const clem::AverageDiagnostics diag = clem::AverageDown(*ClemPC, level);

  if (verbose != 0) {
    amrex::Print() << "CLEM: level " << level << " average down from level " << level + 1 << ": "
                   << diag.n_covered << " covered cells, " << diag.n_rebuilt << " rebuilt"
                   << std::endl;
  }
}
