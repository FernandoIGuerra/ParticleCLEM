/**
 * post_regrid.cpp: the CLEM side of an AMR regrid inside PeleC.
 *
 *   PeleC::SetNewCellsFromCoarse  called when a level is rebuilt on new
 *                                 grids: its new cells take the LES state of
 *                                 their coarse parent, piecewise constant
 *   PeleC::PostRegridClem         called by PeleC::post_regrid of the base
 *                                 level: the lines follow the new grids
 *                                 (grid_adapter.H)
 *
 * The two go together: a new fine cell gets the clone of the line of its
 * parent, so its LES state must be the one of its parent too.
 */
#include <AMReX_MultiFab.H>

#include "PeleC.H"
#include "IndexDefines.H"

#include "grid_adapter.H"
#include "setup.H"
#include "splicer.H"

void
PeleC::SetNewCellsFromCoarse(
  const amrex::BoxArray* old_grids, const amrex::DistributionMapping* old_dmap)
{
  BL_PROFILE("PeleC::SetNewCellsFromCoarse()");
  AMREX_ALWAYS_ASSERT(level > 0);

  amrex::MultiFab& S_new         = get_new_data(State_Type);
  const amrex::IntVect ref_ratio = parent->refRatio(level - 1);

  // The state of the parents, on the rank of their children
  amrex::BoxArray parent_grids = grids;
  parent_grids.coarsen(ref_ratio);
  amrex::MultiFab parents(parent_grids, dmap, NVAR, 0);
  parents.ParallelCopy(
    getLevel(level - 1).get_new_data(State_Type), 0, 0, NVAR, 0, 0,
    parent->Geom(level - 1).periodicity());

  // 1 in the cells that were not on this level before the regrid
  amrex::MultiFab is_new(grids, dmap, 1, 0);
  is_new.setVal(1.0);
  if (old_grids != nullptr) {
    amrex::MultiFab kept(*old_grids, *old_dmap, 1, 0);
    kept.setVal(0.0);
    is_new.ParallelCopy(kept, 0, 0, 1);
  }

  for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const auto& state    = S_new.array(mfi);
    const auto& src      = parents.const_array(mfi);
    const auto& mark     = is_new.const_array(mfi);

    amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
      if (mark(i, j, k) < 0.5) {
        return;
      }
      const amrex::IntVect iv_crse = amrex::coarsen(amrex::IntVect(AMREX_D_DECL(i, j, k)), ref_ratio);
      for (int n = 0; n < NVAR; ++n) {
        state(i, j, k, n) = src(iv_crse, n);
      }
    });
  }
}

void
PeleC::PostRegridClem(const int lbase)
{
  BL_PROFILE("PeleC::PostRegridClem()");
  // The first regrids happen while the grids are built, before the lines exist
  if (ClemPC == nullptr || !ClemPC->IsDefined()) {
    return;
  }

  const clem::AdaptDiagnostics diag = clem::AdaptToGrids(*ClemPC, lbase);

  if (verbose != 0) {
    amrex::Print() << "CLEM: regrid from level " << lbase << ": " << diag.n_dropped
                   << " elements dropped, " << diag.n_refined << " new fine cells cloned"
                   << std::endl;
  }
  if (clem::setup::check_mass_flux) {
    for (int lev = lbase; lev <= parent->finestLevel(); ++lev) {
      amrex::Print() << "CLEM: level " << lev
                     << " line mass defect after the regrid, max |m_line / (rho V) - 1| = "
                     << clem::LineMassDefect(
                          *ClemPC, lev, getLevel(lev).get_new_data(State_Type), URHO)
                     << std::endl;
    }
  }
}
