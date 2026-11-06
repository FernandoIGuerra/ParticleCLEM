/**
 * restart_clem.cpp: checkpoint and restart of the CLEM lines.
 *
 *   ClemContainer::Checkpoint  every element of every level into the
 *                              checkpoint directory, with CLEMHeader
 *   ClemContainer::Restart     the elements back, on the restored grids;
 *                              lines of another n_lem remapped to n_lem
 *   ClemContainer::PadLines    massless tail elements, so that a line of
 *                              fewer than n_lem elements can be regridded
 *
 * What is saved and what is not: the elements are the whole state of the
 * model. The checkpoint is written at the end of a coarse step, after every
 * coarse/fine sync, so the face mass fluxes and transfers are consumed, the
 * halo is refilled at the next step, nu_T is recomputed by the LES closure
 * and the CVODE reactor is a solver context: none of them is saved. The
 * random state of the stirring is not saved either: a restarted run draws
 * other eddies than an uninterrupted one, the same distribution.
 */
#include <fstream>
#include <sstream>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include "averager.H"
#include "particle_container.H"
#include "setup.H"

namespace clem {

namespace {

constexpr int checkpoint_version = 1;
const std::string header_file    = "CLEMHeader";

// the elements of one level, as a single-level particle set
using SingleLevelContainer = amrex::ParticleContainer<RealData::ncomps, IntData::ncomps, 0, 0>;

// chk_dir/clem_particles_level_<lev>
std::string
LevelDir(const int lev)
{
  return "clem_particles_level_" + std::to_string(lev);
}

} // namespace

void
ClemContainer::Checkpoint(const std::string& chk_dir) const
{
  BL_PROFILE("clem::ClemContainer::Checkpoint()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  amrex::Vector<std::string> real_names;
  amrex::Vector<std::string> int_names;
  ParticleComponentNames(real_names, int_names);

  // one single-level particle set per level: an AMR particle container files
  // every particle at the finest level that holds its position, which would
  // move the lines of the covered cells onto the level above at the restart.
  // A single-level container cannot: each line stays on its level
  for (int lev = 0; lev <= m_pc->finestLevel(); ++lev) {
    SingleLevelContainer one(m_pc->Geom(lev), m_pc->ParticleDistributionMap(lev),
                             m_pc->ParticleBoxArray(lev));
    for (const auto& [key, tile] : m_pc->GetParticles(lev)) {
      auto& dst = one.DefineAndReturnParticleTile(0, key.first, key.second);
      for (const auto& p : tile.GetArrayOfStructs()) {
        dst.push_back(p);
      }
    }
    one.Checkpoint(chk_dir, LevelDir(lev), real_names, int_names);
  }

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
    std::ofstream os(chk_dir + "/" + header_file);
    os << "CLEM checkpoint version " << checkpoint_version << "\n";
    os << "n_lem " << setup::n_lem << "\n";
    os << "species " << spec_names.size();
    for (const auto& name : spec_names) {
      os << " " << name;
    }
    os << "\n";
    os << "finest_level " << m_pc->finestLevel() << "\n";
  }
}

bool
ClemContainer::Restart(const std::string& chk_dir)
{
  BL_PROFILE("clem::ClemContainer::Restart()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  if (!amrex::FileExists(chk_dir + "/" + header_file)) {
    return false;
  }

  // the header: the line size it was written with, and the species, which
  // must be those of this build
  amrex::Vector<char> buffer;
  amrex::ParallelDescriptor::ReadAndBcastFile(chk_dir + "/" + header_file, buffer);
  std::istringstream is(std::string(buffer.dataPtr()));
  std::string word;
  int version = 0;
  int n_saved = 0;
  int n_spec  = 0;
  is >> word >> word >> word >> version;      // CLEM checkpoint version <v>
  is >> word >> n_saved;                      // n_lem <n>
  is >> word >> n_spec;                       // species <n> <names>
  amrex::Vector<std::string> saved_names(n_spec);
  for (auto& name : saved_names) {
    is >> name;
  }
  int saved_finest = -1;
  is >> word >> saved_finest;                 // finest_level <l>
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(version == checkpoint_version && n_saved > 0,
                                   "CLEM restart: unreadable CLEMHeader");
  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(saved_names == spec_names,
                                   "CLEM restart: the checkpoint has other species than this build");

  // a restart has at most the levels of its checkpoint (AMReX restores the
  // saved hierarchy; a lower amr.max_level drops the finest ones)
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_pc->finestLevel() <= saved_finest,
                                   "CLEM restart: the hierarchy has more levels than the checkpoint");

  // each level read into a single-level set on the restored grids of that
  // level, whatever the rank count; its restart files every element in the
  // right box of those grids, then the boxes are moved into the AMR container
  m_pc->resizeData();
  for (int lev = 0; lev <= m_pc->finestLevel(); ++lev) {
    SingleLevelContainer one(m_pc->Geom(lev), m_pc->ParticleDistributionMap(lev),
                             m_pc->ParticleBoxArray(lev));
    one.Restart(chk_dir, LevelDir(lev));
    for (auto& [key, tile] : one.GetParticles(0)) {
      auto& dst = m_pc->DefineAndReturnParticleTile(lev, key.first, key.second);
      for (const auto& p : tile.GetArrayOfStructs()) {
        dst.push_back(p);
      }
    }
  }
  InvalidateCellMaps();

  if (n_saved != setup::n_lem) {
    amrex::Print() << "CLEM restart: lines of " << n_saved << " elements remapped to " << setup::n_lem
                   << std::endl;
    const setup::RegridderType regridder;
    for (int lev = 0; lev <= m_pc->finestLevel(); ++lev) {
      if (n_saved < setup::n_lem) {
        PadLines(lev, n_saved);
      }
      regridder.Regrid(*this, lev);
    }
    // element k of a covered cell is the average of element k of its children
    AverageDownAll(*this);
  }
  RebuildCellMaps();
  return true;
}

void
ClemContainer::PadLines(const int lev, const int n_old)
{
  if (!IsCellMapValid(lev)) {
    RebuildCellMap(lev);
  }
  const amrex::Geometry& geom = m_pc->Geom(lev);
  const auto prob_lo          = geom.ProbLoArray();
  const auto cell_size        = geom.CellSizeArray();

  std::vector<ParticleType> fresh;
  for (amrex::MFIter mfi = m_pc->MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx  = mfi.validbox();
    auto& tile            = m_pc->DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
    const auto* particles = tile.GetArrayOfStructs().data();

    fresh.clear();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& in_cell = ParticlesInCell(lev, mfi.index(), iv);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(in_cell.size()) == n_old,
                                       "CLEM restart: a line does not have the n_lem of its checkpoint");
      // the tail element: the template of the massless ones
      int tail = in_cell[0];
      for (const int n : in_cell) {
        if (particles[n].idata(IntData::lem_index) > particles[tail].idata(IntData::lem_index)) {
          tail = n;
        }
      }
      for (int k = n_old; k < setup::n_lem; ++k) {
        ParticleType p = particles[tail];
        p.id()         = ParticleType::NextID();
        p.cpu()        = amrex::ParallelDescriptor::MyProc();
        p.rdata(RealData::mass) = 0.0;
        p.rdata(RealData::vol)  = 0.0;
        PlaceOnLine(p, iv, k, setup::n_lem, prob_lo, cell_size);
        fresh.push_back(p);
      }
    }
    // after the loop: push_back may move the storage of the tile
    for (const auto& p : fresh) {
      tile.push_back(p);
    }
  }
  m_cell_map_valid[lev] = 0;
}

} // namespace clem
