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
    pp.query("phi", parm.phi);
    pp.query("p0", parm.p0);
    pp.query("T0", parm.T0);
    pp.query("T_ign", parm.T_ign);
    pp.query("r_ign", parm.r_ign);
  }

  // H2/air: phi H2 + 0.5 (O2 + 3.76 N2), as mole fractions
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  const amrex::Real n_h2 = parm.phi;
  const amrex::Real n_o2 = 0.5;
  const amrex::Real n_n2 = 0.5 * 3.76;
  const amrex::Real n_tot = n_h2 + n_o2 + n_n2;
  molefrac[H2_ID] = n_h2 / n_tot;
  molefrac[O2_ID] = n_o2 / n_tot;
  molefrac[N2_ID] = n_n2 / n_tot;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  for (int n = 0; n < NUM_SPECIES; n++) {
    parm.massfrac[n] = massfrac[n];
  }

  // the mixture and the kernel, both at p0
  eos.PYT2RE(parm.p0, massfrac, parm.T0, parm.rho0, parm.e0);
  eos.PYT2RE(parm.p0, massfrac, parm.T_ign, parm.rho_ign, parm.e_ign);

  amrex::Print() << "ClemTulipFlame: H2/air phi = " << parm.phi << " at " << parm.p0 << ", "
                 << parm.T0 << " K; kernel r = " << parm.r_ign << " cm at " << parm.T_ign << " K"
                 << std::endl;
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
