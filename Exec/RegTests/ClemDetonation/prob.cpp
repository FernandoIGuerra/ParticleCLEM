#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  ProbParmDevice& parm = *PeleC::h_prob_parm_device;
  {
    amrex::ParmParse pp("prob");
    pp.query("X_H2", parm.X_H2);
    pp.query("X_O2", parm.X_O2);
    pp.query("X_AR", parm.X_AR);
    pp.query("p0", parm.p0);
    pp.query("T0", parm.T0);
    pp.query("p_driver", parm.p_driver);
    pp.query("T_driver", parm.T_driver);
    pp.query("x_driver", parm.x_driver);
  }

  // the mixture: mole fractions -> mass fractions
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  const amrex::Real sum = parm.X_H2 + parm.X_O2 + parm.X_AR;
  molefrac[H2_ID] = parm.X_H2 / sum;
  molefrac[O2_ID] = parm.X_O2 / sum;
  molefrac[AR_ID] = parm.X_AR / sum;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  for (int n = 0; n < NUM_SPECIES; n++) {
    parm.massfrac[n] = massfrac[n];
  }

  // the two states from (p, T, Y)
  eos.PYT2RE(parm.p0, massfrac, parm.T0, parm.rho0, parm.e0);
  eos.PYT2RE(parm.p_driver, massfrac, parm.T_driver, parm.rho_driver, parm.e_driver);

  amrex::Print() << "ClemDetonation: mixture H2:O2:AR = " << parm.X_H2 << ":" << parm.X_O2 << ":"
                 << parm.X_AR << ", rho0 = " << parm.rho0 << " g/cm3 at p0 = " << parm.p0
                 << ", T0 = " << parm.T0 << "; driver x < " << parm.x_driver << " cm at "
                 << parm.p_driver << ", " << parm.T_driver << " K" << std::endl;
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
