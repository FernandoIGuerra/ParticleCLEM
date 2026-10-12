/**
 * initialize_clem.cpp: the start of the CLEM model inside PeleC.
 *
 *   PeleC::InitializeClem  build the particle container on the AMR hierarchy
 *                          and seed the lines of every level from the LES
 */
#include "PeleC.H"
#include "IndexDefines.H"

#include "averager.H"
#include "particle_container.H"

void
PeleC::InitializeClem()
{
  BL_PROFILE("PeleC::InitializeClem()");
  AMREX_ALWAYS_ASSERT(level == 0 && ClemPC == nullptr);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    do_mol && mol_iters == 1,
    "CLEM needs the MOL integrator with one iteration: pelec.do_mol = 1, pelec.mol_iters = 1");

  ClemPC = std::make_unique<clem::ClemContainer>();
  ClemPC->Define(parent);

  // called from post_init, once the grids of every level are final and the
  // LES state is averaged down: each level is seeded from its own state
  for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
    const amrex::MultiFab& les = getLevel(lev).get_new_data(State_Type);
    ClemPC->InitializeFromLES(lev, les, URHO, UTEMP, UEINT, UFS);

    if (verbose != 0) {
      amrex::Print() << "CLEM: level " << lev << " seeded with "
                     << ClemPC->Particles().NumberOfParticlesAtLevel(lev) << " particles in "
                     << les.boxArray().numPts() << " cells" << std::endl;
    }
  }

  // every covered cell starts as the average of its children, as it will be
  // after every sync
  clem::AverageDownAll(*ClemPC);
}
