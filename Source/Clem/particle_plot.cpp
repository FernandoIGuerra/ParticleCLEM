/**
 * particle_plot.cpp: the plotfiles of the LEM lines (ClemContainer).
 *
 *   WritePlotFile      every particle of every level
 *   WriteLeafPlotFile  the particles of the leaves only
 */
#include <cctype>
#include <fstream>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include "PelePhysics.H"

#include "particle_container.H"

namespace clem {

namespace {
// the names of the particle components, in the order of RealData / IntData
void
ComponentNames(amrex::Vector<std::string>& real_names, amrex::Vector<std::string>& int_names)
{
  real_names.resize(RealData::ncomps);
  real_names[RealData::rho]   = "rho";
  real_names[RealData::mass]  = "mass";
  real_names[RealData::vol]   = "vol";
  real_names[RealData::T]     = "T";
  real_names[RealData::eint]  = "eint";
  real_names[RealData::press] = "press";
  real_names[RealData::phi]   = "phi";

  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  
  for (int n = 0; n < NUM_SPECIES; ++n) {
    real_names[RealData::Y0 + n] = "Y(" + spec_names[n] + ")";
  }

  int_names.resize(IntData::ncomps);
  int_names[IntData::lem_index]       = "lem_index";
  int_names[IntData::entry_side]      = "entry_side";
}

// the particle directory of a LES plotfile: next to it, the tag between the
// root (amr.plot_file) and the step, plt00048 -> plt<tag>00048
std::string
ParticlePlotFileName(std::string grid_plotfile, const std::string& tag)
{
  while (!grid_plotfile.empty() && grid_plotfile.back() == '/') {
    grid_plotfile.pop_back();
  }
  // the step: the trailing digits
  auto digits = grid_plotfile.size();
  while (digits > 0 && std::isdigit(static_cast<unsigned char>(grid_plotfile[digits - 1])) != 0) {
    --digits;
  }
  return grid_plotfile.substr(0, digits) + tag + grid_plotfile.substr(digits);
}

// write the particles of pc into particle_dir/particles, and copy the header
// of the LES plotfile: it makes particle_dir a plotfile ParaView accepts
void
WriteParticles(const ClemParticleContainer& pc, const std::string& grid_plotfile,
               const std::string& particle_dir)
{
  amrex::Vector<std::string> real_names;
  amrex::Vector<std::string> int_names;

  ComponentNames(real_names, int_names);
  pc.WritePlotFile(particle_dir, "particles", real_names, int_names);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ifstream src(grid_plotfile + "/Header", std::ios::binary);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(src.good(), "CLEM: the LES plotfile has no Header yet");
    std::ofstream dst(particle_dir + "/Header", std::ios::binary | std::ios::trunc);
    dst << src.rdbuf();
  }
}
} // namespace

void
ClemContainer::WritePlotFile(const std::string& grid_plotfile) const
{
  BL_PROFILE("clem::ClemContainer::WritePlotFile()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  WriteParticles(*m_pc, grid_plotfile, ParticlePlotFileName(grid_plotfile, "_particle_"));
}

void
ClemContainer::WriteLeafPlotFile(const std::string& grid_plotfile)
{
  BL_PROFILE("clem::ClemContainer::WriteLeafPlotFile()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  FlipCoveredParticles(); // hide
  WriteParticles(*m_pc, grid_plotfile, ParticlePlotFileName(grid_plotfile, "_particle_leave_"));
  FlipCoveredParticles(); // restore
}

void
ClemContainer::FlipCoveredParticles()
{
  // the finest level covers nothing
  for (int lev = 0; lev < m_pc->finestLevel(); ++lev) {
    // 1 where the cell is covered by level lev + 1, 0 where it is a leaf
    const amrex::iMultiFab covered = amrex::makeFineMask(
      m_pc->ParticleBoxArray(lev), m_pc->ParticleDistributionMap(lev),
      m_pc->ParticleBoxArray(lev + 1), m_pc->GetParGDB()->refRatio(lev), 0, 1);

    for (ClemParIter pti(*m_pc, lev); pti.isValid(); ++pti) {

      auto& aos           = pti.GetArrayOfStructs();
      const auto& mask    = covered.const_array(pti.index());

      for (int n = 0; n < pti.numParticles(); ++n) {
        auto& p = aos[n];
        if (mask(m_pc->Index(p, lev)) == 1) {
          p.id() = -p.id();
        }
      }
    }
  }
}

} // namespace clem
