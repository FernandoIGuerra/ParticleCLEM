#include <AMReX_ParallelDescriptor.H>

#include "PelePhysics.H"

#include "particle_container.H"
#include "setup.H"

namespace clem {

void
ClemContainer::Define(amrex::AmrCore* amr)
{
  AMREX_ALWAYS_ASSERT(m_pc == nullptr);
  m_pc = std::make_unique<ClemParticleContainer>(amr);
  InvalidateCellMaps();
}

void
ClemContainer::InitializeFromLES(const int lev, const amrex::MultiFab& les, const int rho_comp,
                                 const int temp_comp, const int eint_comp, const int spec_comp)
{
  BL_PROFILE("clem::ClemContainer::InitializeFromLES()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    les.boxArray() == m_pc->ParticleBoxArray(lev) && les.DistributionMap() == m_pc->ParticleDistributionMap(lev),
    "InitializeFromLES: LES state must live on the particle grids of this level");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_pc->NumberOfParticlesAtLevel(lev) == 0,
                                   "InitializeFromLES: level already holds particles");

  const amrex::Geometry& geom = m_pc->Geom(lev);
  const auto prob_lo          = geom.ProbLoArray();
  const auto cell_size        = geom.CellSizeArray();
  const amrex::Real cell_vol  = AMREX_D_TERM(cell_size[0], *cell_size[1], *cell_size[2]);
  
  // equal volume share per element, so sum_k m_k = rho V_cell
  const amrex::Real vol = cell_vol / static_cast<amrex::Real>(setup::n_lem);
  // element k at fraction (k+1)/(n+1) of the cell along x, strictly interior,
  // so AMReX bins it in its own cell
  const amrex::Real dx_inner = cell_size[0] / static_cast<amrex::Real>(setup::n_lem + 1);

  auto eos = pele::physics::PhysicsType::eos();

  for (amrex::MFIter mfi = m_pc->MakeMFIter(lev, false); mfi.isValid(); ++mfi) {

    const amrex::Box& bx  = mfi.validbox();
    const auto& state     = les.const_array(mfi);
    auto& tile            = m_pc->DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
    tile.GetArrayOfStructs()().reserve(bx.numPts() * setup::n_lem);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {

      const amrex::Real rho = state(iv, rho_comp);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rho > 0.0, "InitializeFromLES: rho <= 0");

      // the LES state stores rho e and rho Y_k: convert to specific
      const amrex::Real temp  = state(iv, temp_comp);
      const amrex::Real eint  = state(iv, eint_comp) / rho;
      
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = state(iv, spec_comp + n) / rho;
      }
      // element pressure = EOS pressure of its own state (the cell pressure)
      amrex::Real press = 0.0;
      eos.RTY2P(rho, temp, massfrac, press);

      for (int k = 0; k < setup::n_lem; ++k) {
        ParticleType p;
        p.id()  = ParticleType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();
        AMREX_D_TERM(
                    p.pos(0) = prob_lo[0] + iv[0] * cell_size[0] + static_cast<amrex::Real>(k + 1) * dx_inner;
                  , p.pos(1) = prob_lo[1] + (iv[1] + 0.5) * cell_size[1];
                  , p.pos(2) = prob_lo[2] + (iv[2] + 0.5) * cell_size[2];)
        p.rdata(RealData::rho)   = rho;
        p.rdata(RealData::mass)  = rho * vol;
        p.rdata(RealData::vol)   = vol;
        p.rdata(RealData::T)     = temp;
        p.rdata(RealData::eint)  = eint;
        p.rdata(RealData::press) = press;
        p.rdata(RealData::phi)   = 0.0;
        for (int n = 0; n < NUM_SPECIES; ++n) {
          p.rdata(RealData::Y0 + n) = massfrac[n];
        }
        p.idata(IntData::lem_index)  = k;
        p.idata(IntData::entry_side) = 0;
        tile.push_back(p);
      }
    }
  }
  // every element is interior to its own cell on its own grid, so no
  // Redistribute is needed; the map of this level is rebuilt directly
  RebuildCellMap(lev);
}

void
ClemContainer::GetFilteredScalars(const int lev, amrex::MultiFab& filtered) const
{
  BL_PROFILE("clem::ClemContainer::GetFilteredScalars()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    filtered.boxArray() == m_pc->ParticleBoxArray(lev) && filtered.DistributionMap() == m_pc->ParticleDistributionMap(lev),
    "GetFilteredScalars: filtered must live on the particle grids of this level");

  AMREX_ALWAYS_ASSERT(filtered.nComp() >= RealData::ncomps);

  // sum m, sum v and sum m q per cell, straight from the particles, so the
  // result does not depend on the cell map being current
  amrex::MultiFab sums(filtered.boxArray(), filtered.DistributionMap(), RealData::ncomps, 0);
  sums.setVal(0.0);

  for (ClemParConstIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    const auto& aos = pti.GetArrayOfStructs();
    const auto& sum = sums.array(pti.index());

    for (int n = 0; n < pti.numParticles(); ++n) {
      const auto& p             = aos[n];
      const amrex::IntVect iv   = m_pc->Index(p, lev);
      const amrex::Real mass    = p.rdata(RealData::mass);
      sum(iv, RealData::mass)   += mass;
      sum(iv, RealData::vol)    += p.rdata(RealData::vol);

      for (int comp = RealData::T; comp < RealData::ncomps; ++comp) {
        sum(iv, comp) += mass * p.rdata(comp);
      }
    }
  }

  for (amrex::MFIter mfi(filtered); mfi.isValid(); ++mfi) {
    
    const amrex::Box& bx = mfi.validbox();
    const auto& sum      = sums.const_array(mfi);
    const auto& out      = filtered.array(mfi);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const amrex::Real mass = sum(iv, RealData::mass);
      const amrex::Real vol  = sum(iv, RealData::vol);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mass > 0.0 && vol > 0.0, "GetFilteredScalars: a cell holds no mass");
      
      out(iv, RealData::mass) = mass;
      out(iv, RealData::vol)  = vol;
      out(iv, RealData::rho)  = mass / vol;

      for (int comp = RealData::T; comp < RealData::ncomps; ++comp) {
        out(iv, comp) = sum(iv, comp) / mass;
      }
      
    }
  }
}

void
ClemContainer::RememberGrids()
{
  const int n_levels = m_pc->finestLevel() + 1;
  m_grids.resize(n_levels);
  m_dmaps.resize(n_levels);
  for (int lev = 0; lev < n_levels; ++lev) {
    m_grids[lev] = m_pc->ParticleBoxArray(lev);
    m_dmaps[lev] = m_pc->ParticleDistributionMap(lev);
  }
}

bool
ClemContainer::GridsChanged(const int lev) const
{
  return lev >= static_cast<int>(m_grids.size()) || m_grids[lev] != m_pc->ParticleBoxArray(lev) ||
         m_dmaps[lev] != m_pc->ParticleDistributionMap(lev);
}

void
ClemContainer::InvalidateCellMaps()
{
  for (auto&& valid : m_cell_map_valid) {
    valid = 0;
  }
}

/**
 * Rebuils particle indexing: LES ownership <--> particle
 */
void
ClemContainer::RebuildCellMap(const int lev)
{
  BL_PROFILE("clem::ClemContainer::RebuildCellMap()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  if (static_cast<int>(m_cell_vectors.size()) < lev + 1) {
    m_cell_vectors.resize(lev + 1);
    m_cell_map_valid.resize(lev + 1, 0);
  }

  // reset the lists on the valid boxes of this level; grids removed by a
  // regrid are dropped together with the old map
  m_cell_vectors[lev].clear();
  for (amrex::MFIter mfi = m_pc->MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    auto& fab = m_cell_vectors[lev][mfi.index()];
    fab.resize(mfi.validbox());
    fab.setVal(std::vector<int>{});
  }

  // bin each particle into its cell, as an index into the grid AoS; untiled
  // layout, so that index is grid-global
  for (ClemParIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    AMREX_ALWAYS_ASSERT(pti.LocalTileIndex() == 0);
    const auto& aos = pti.GetArrayOfStructs();
    auto& fab       = m_cell_vectors[lev][pti.index()];
    for (int n = 0; n < pti.numParticles(); ++n) {
      fab(m_pc->Index(aos[n], lev)).push_back(n);
    }
  }
  m_cell_map_valid[lev] = 1;
}

void
ClemContainer::RebuildCellMaps()
{
  for (int lev = 0; lev <= m_pc->finestLevel(); ++lev) {
    RebuildCellMap(lev);
  }
}

// ---------------------------------------------------------------------------------------------------------
const std::vector<int>&
ClemContainer::ParticlesInCell(const int lev, const int grid_id, const amrex::IntVect& iv) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(IsCellMapValid(lev), "ClemContainer: cell map used while stale");
  return m_cell_vectors[lev].at(grid_id)(iv);
}

/**
 * Redistribution of particles, between ranks
 */
void
ClemContainer::Redistribute(const int lev)
{
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);
  m_pc->Redistribute(lev, lev);
  if (lev < static_cast<int>(m_cell_map_valid.size())) {
    m_cell_map_valid[lev] = 0;
  }
}

void
ClemContainer::Redistribute()
{
  for (int lev = 0; lev <= m_pc->finestLevel(); ++lev) {
    Redistribute(lev);
  }
}

} // namespace clem
