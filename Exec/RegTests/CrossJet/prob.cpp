#include "prob.H"

namespace {

// Index of a species by name, so the setup is not tied to one mechanism.
int
species_index(const amrex::Vector<std::string>& names, const std::string& s)
{
  for (int n = 0; n < names.size(); n++) {
    if (names[n] == s) {
      return n;
    }
  }
  amrex::Abort("CrossJet: species " + s + " is not in the chemistry model");
  return -1;
}

} // namespace

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
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  auto& pp_d = *PeleC::h_prob_parm_device;

  amrex::ParmParse pp("prob");
  pp.query("p_in", pp_d.p_in);
  pp.query("T_in", pp_d.T_in);
  pp.query("u_in", pp_d.u_in);
  pp.query("YO2_in", pp_d.YO2_in);
  pp.query("YN2_in", pp_d.YN2_in);
  pp.query("p_jet", pp_d.p_jet);
  pp.query("T_jet", pp_d.T_jet);
  pp.query("v_jet", pp_d.v_jet);
  pp.query("jet_x", pp_d.jet_x);
  pp.query("jet_z", pp_d.jet_z);
  pp.query("jet_d", pp_d.jet_d);
  pp.query("jet_start", pp_d.jet_start);
  pp.query("jet_ramp", pp_d.jet_ramp);

  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    spec_names);
  const int iH2 = species_index(spec_names, "H2");
  const int iO2 = species_index(spec_names, "O2");
  const int iN2 = species_index(spec_names, "N2");

  for (int n = 0; n < NUM_SPECIES; n++) {
    pp_d.Y_in[n] = 0.0;
    pp_d.Y_jet[n] = 0.0;
  }
  pp_d.Y_in[iO2] = pp_d.YO2_in;
  pp_d.Y_in[iN2] = pp_d.YN2_in;
  pp_d.Y_jet[iH2] = 1.0;

  amrex::Real Y_in[NUM_SPECIES] = {0.0};
  amrex::Real Y_jet[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y_in[n] = pp_d.Y_in[n];
    Y_jet[n] = pp_d.Y_jet[n];
  }

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(pp_d.p_in, Y_in, pp_d.T_in, pp_d.rho_in, pp_d.e_in);
  eos.PYT2RE(pp_d.p_jet, Y_jet, pp_d.T_jet, pp_d.rho_jet, pp_d.e_jet);

  // Report the derived state so it can be checked against Table 1 of
  // Jin et al. (2024) before committing to a long run.
  amrex::Real cs_in = 0.0;
  amrex::Real cs_jet = 0.0;
  eos.RTY2Cs(pp_d.rho_in, pp_d.T_in, Y_in, cs_in);
  eos.RTY2Cs(pp_d.rho_jet, pp_d.T_jet, Y_jet, cs_jet);

  amrex::Real mw[NUM_SPECIES] = {0.0};
  eos.molecular_weight(mw);
  // 2 H2 + O2 -> 2 H2O
  const amrex::Real fs = 2.0 * mw[iH2] / mw[iO2];

#if AMREX_SPACEDIM == 3
  const amrex::Real inlet_area =
    (probhi[1] - problo[1]) * (probhi[2] - problo[2]);
  const amrex::Real jet_area =
    0.25 * M_PI * pp_d.jet_d * pp_d.jet_d; // cm^2
  const std::string per = "g/s";
#else
  // Per unit spanwise depth in 2D.
  const amrex::Real inlet_area = (probhi[1] - problo[1]);
  const amrex::Real jet_area = pp_d.jet_d;
  const std::string per = "g/s/cm";
#endif

  const amrex::Real mdot_air = pp_d.rho_in * pp_d.u_in * inlet_area;
  const amrex::Real mdot_jet = pp_d.rho_jet * pp_d.v_jet * jet_area;
  const amrex::Real phi = (mdot_jet / (pp_d.YO2_in * mdot_air)) / fs;
  const amrex::Real J = (pp_d.rho_jet * pp_d.v_jet * pp_d.v_jet) /
                        (pp_d.rho_in * pp_d.u_in * pp_d.u_in);

  amrex::Print() << "\n=== CrossJet: H2 jet in supersonic crossflow ===\n"
                 << "  crossflow  p       = " << pp_d.p_in / 1.01325e6
                 << " atm\n"
                 << "  crossflow  T       = " << pp_d.T_in << " K\n"
                 << "  crossflow  u       = " << pp_d.u_in * 1.0e-2 << " m/s\n"
                 << "  crossflow  rho     = " << pp_d.rho_in << " g/cm^3\n"
                 << "  crossflow  Mach    = " << pp_d.u_in / cs_in << "\n"
                 << "  crossflow  mdot    = " << mdot_air << " " << per << "\n"
                 << "  jet        p       = " << pp_d.p_jet / 1.01325e6
                 << " atm\n"
                 << "  jet        T       = " << pp_d.T_jet << " K\n"
                 << "  jet        v       = " << pp_d.v_jet * 1.0e-2 << " m/s\n"
                 << "  jet        rho     = " << pp_d.rho_jet << " g/cm^3\n"
                 << "  jet        Mach    = " << pp_d.v_jet / cs_jet << "\n"
                 << "  jet        mdot    = " << mdot_jet << " " << per << "\n"
                 << "  global equiv ratio = " << phi << "\n"
                 << "  momentum ratio J   = " << J << "\n"
                 << "  injector centre    = (" << pp_d.jet_x << ", "
                 << problo[1] << ", " << pp_d.jet_z << ") cm, d = "
                 << pp_d.jet_d << " cm\n"
                 << "================================================\n\n";
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
