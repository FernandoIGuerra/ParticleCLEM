/**
 * initialize_clem.cpp: the start of the CLEM model inside PeleC.
 *
 *   PeleC::InitializeClem  build the particle container on the AMR hierarchy
 *                          and seed the lines of every level from the LES,
 *                          or read them back from a checkpoint at a restart
 */
#include "PeleC.H"
#include "IndexDefines.H"

#include "averager.H"
#include "particle_container.H"
#include "setup.H"
#include "splicer.H"

void
PeleC::InitializeClem(const std::string& restart_dir)
{
  BL_PROFILE("PeleC::InitializeClem()");
  AMREX_ALWAYS_ASSERT(level == 0 && ClemPC == nullptr);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    do_mol && mol_iters == 1,
    "CLEM needs the MOL integrator with one iteration: pelec.do_mol = 1, pelec.mol_iters = 1");

  ClemPC = std::make_unique<clem::ClemContainer>();
  ClemPC->Define(parent);

  // at a restart (post_restart, every level restored) the lines come back
  // from the checkpoint; a checkpoint without them (pure PeleC, or written
  // before CLEM had a restart) seeds them from its LES state, as a new run
  bool restored = false;
  if (!restart_dir.empty()) {
    restored = ClemPC->Restart(restart_dir);
    amrex::Print() << (restored ? "CLEM: lines read from " : "CLEM: no lines in ") << restart_dir
                   << (restored ? "" : ", seeded from its LES state") << std::endl;
    for (int lev = 0; restored && lev <= parent->finestLevel(); ++lev) {
      const amrex::Long n_read = ClemPC->Particles().NumberOfParticlesAtLevel(lev);
      const amrex::Long n_want = getLevel(lev).boxArray().numPts() * clem::setup::n_lem;
      amrex::Print() << "CLEM: level " << lev << " restored with " << n_read << " particles, n_lem x cells = "
                     << n_want << std::endl;
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_read == n_want,
                                       "CLEM restart: a level does not hold n_lem elements per cell");
      if (clem::setup::check_mass_flux) {
        amrex::Print() << "CLEM: level " << lev
                       << " line mass defect after the restart, max |m_line / (rho V) - 1| = "
                       << clem::LineMassDefect(*ClemPC, lev, getLevel(lev).get_new_data(State_Type), URHO)
                       << std::endl;
      }
    }
  }

  // called from post_init, once the grids of every level are final and the
  // LES state is averaged down: each level is seeded from its own state
  for (int lev = 0; !restored && lev <= parent->finestLevel(); ++lev) {
    const amrex::MultiFab& les = getLevel(lev).get_new_data(State_Type);
    ClemPC->InitializeFromLES(lev, les, URHO, UTEMP, UEINT, UFS);

    if (verbose != 0) {
      amrex::Print() << "CLEM: level " << lev << " seeded with "
                     << ClemPC->Particles().NumberOfParticlesAtLevel(lev) << " particles in "
                     << les.boxArray().numPts() << " cells" << std::endl;
    }
  }

  // every covered cell starts as the average of its children, as it will be
  // after every sync (a checkpoint is written after the sync: already so)
  if (!restored) {
    clem::AverageDownAll(*ClemPC);
  }

  // the grids the lines are filed under, to tell when an AMR regrid moves them
  ClemPC->RememberGrids();

  // the reactor of the elements: CLEM works only with CVODE
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(chem_integrator == "ReactorCvode",
                                   "CLEM needs pelec.chem_integrator = ReactorCvode");
  ClemPC->InitializeReactor();

  // its typical values follow PeleC's (use_typ_vals_chem): set now, refreshed
  // whenever PeleC refreshes its own
  if (use_typical_vals_chem) {
    set_typical_values_chem();
  }
}
