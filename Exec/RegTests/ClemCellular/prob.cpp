#include <fstream>
#include <sstream>

#include "prob.H"

void
pc_prob_close()
{
}

namespace {

// the ZND table: a comment line, the number of points, then per point the
// distance behind the shock and the ZND_NCOL columns
void
ReadProfile(const std::string& name)
{
  std::ifstream in(name);
  if (!in.good()) {
    amrex::Abort("ClemCellular: cannot open the ZND table " + name);
  }
  std::string line;
  std::getline(in, line);
  std::getline(in, line);
  const int n = std::stoi(line);

  auto& host = *PeleC::prob_parm_host;
  host.h_znd_x.resize(n);
  host.h_znd_q.resize(static_cast<long>(n) * ZND_NCOL);
  host.znd_x.resize(n);
  host.znd_q.resize(static_cast<long>(n) * ZND_NCOL);
  for (int i = 0; i < n; i++) {
    in >> host.h_znd_x[i];
    for (int c = 0; c < ZND_NCOL; c++) {
      in >> host.h_znd_q[ZND_NCOL * i + c];
    }
  }
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(in.good(), "ClemCellular: the ZND table is incomplete");

  amrex::Gpu::copy(amrex::Gpu::hostToDevice, host.h_znd_x.begin(), host.h_znd_x.end(), host.znd_x.begin());
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, host.h_znd_q.begin(), host.h_znd_q.end(), host.znd_q.begin());
  PeleC::h_prob_parm_device->znd_n = n;
  PeleC::h_prob_parm_device->d_znd_x = host.znd_x.data();
  PeleC::h_prob_parm_device->d_znd_q = host.znd_q.data();
}

} // namespace

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
    pp.query("p0", parm.p0);
    pp.query("T0", parm.T0);
    pp.query("X_H2", parm.X_H2);
    pp.query("X_O2", parm.X_O2);
    pp.query("X_AR", parm.X_AR);
    pp.query("D_window", parm.D_window);
    pp.query("x_shock", parm.x_shock);
    pp.query("alpha", parm.alpha);
    pp.query("hs_p", parm.hs_p);
    pp.query("hs_x", parm.hs_x);
    pp.query("hs_y", parm.hs_y);
    pp.query("hs_r", parm.hs_r);
    pp.query("profile", PeleC::prob_parm_host->profile);
  }
  ReadProfile(PeleC::prob_parm_host->profile);

  // the fresh mixture: mole fractions -> mass fractions, (p, T) -> (rho, e)
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  const amrex::Real sum             = parm.X_H2 + parm.X_O2 + parm.X_AR;

  molefrac[H2_ID] = parm.X_H2 / sum;
  molefrac[O2_ID] = parm.X_O2 / sum;
  molefrac[AR_ID] = parm.X_AR / sum;

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  
  for (int n = 0; n < NUM_SPECIES; n++) {
    parm.massfrac0[n] = massfrac[n];
  }
  eos.PYT2RE(parm.p0, massfrac, parm.T0, parm.rho0, parm.e0);

  amrex::Print() << "ClemCellular: H2:O2:AR = " << parm.X_H2 << ":" << parm.X_O2 << ":" << parm.X_AR
                 << ", p0 = " << parm.p0 << ", T0 = " << parm.T0 << ", rho0 = " << parm.rho0
                 << "; window at " << parm.D_window << " cm/s, shock at x = " << parm.x_shock
                 << ", disturbance alpha = " << parm.alpha << ", " << parm.znd_n
                 << " points of ZND" << std::endl;
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
