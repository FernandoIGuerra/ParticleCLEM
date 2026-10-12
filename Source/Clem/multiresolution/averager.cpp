#include <cmath>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_iMultiFab.H>

#include "PelePhysics.H"

#include "averager.H"
#include "eos_util.H"
#include "setup.H"

namespace clem {

namespace {

using ParticleType = ClemContainer::ParticleType;

// what is summed over the children, per element: m, vol, m e, m T, m phi,
// m Y_k. m T only gives the first guess of the temperature
constexpr int c_mass      = 0;
constexpr int c_vol       = 1;
constexpr int c_mass_eint = 2;
constexpr int c_mass_T    = 3;
constexpr int c_mass_phi  = 4;
constexpr int c_mass_Y    = 5;
constexpr int n_sums      = c_mass_Y + NUM_SPECIES;

// equal mass within a line, as left by the regridding
constexpr amrex::Real equal_mass_tol = 1.0e-9;

// the particle of a cell with lem_index k, or -1
int
FindElement(const std::vector<int>& in_cell, const ParticleType* particles, const int k)
{
  for (const int n : in_cell) {
    if (particles[n].idata(IntData::lem_index) == k) {
      return n;
    }
  }
  return -1;
}

// does the cell hold exactly one particle for each lem_index 0 .. n_lem - 1
bool
HasAllSlots(const std::vector<int>& in_cell, const ParticleType* particles)
{
  if (static_cast<int>(in_cell.size()) != setup::n_lem) {
    return false;
  }
  for (int k = 0; k < setup::n_lem; ++k) {
    if (FindElement(in_cell, particles, k) < 0) {
      return false;
    }
  }
  return true;
}

// every covered cell of level lev gets its n_lem slots, lem_index 0 ..
// n_lem - 1: kept if it already holds them, replaced otherwise. `covered` is
// 1 in the covered cells
void
PrepareCoveredCells(ClemContainer& container, const int lev, const amrex::iMultiFab& covered,
                    AverageDiagnostics& diag)
{
  auto& pc             = container.Particles();
  const auto prob_lo   = pc.Geom(lev).ProbLoArray();
  const auto cell_size = pc.Geom(lev).CellSizeArray();

  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev);
  }

  std::vector<ParticleType> fresh;
  amrex::Long n_invalid = 0;
  for (amrex::MFIter mfi = pc.MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    const int grid_id    = mfi.index();
    auto& tile           = pc.DefineAndReturnParticleTile(lev, grid_id, mfi.LocalTileIndex());
    auto* particles      = tile.GetArrayOfStructs().data();
    const auto& mask     = covered.const_array(mfi);

    fresh.clear();
    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (mask(iv) != 1) {
        continue;
      }
      ++diag.n_covered;
      const auto& in_cell = container.ParticlesInCell(lev, grid_id, iv);
      if (HasAllSlots(in_cell, particles)) {
        continue;
      }
      ++diag.n_rebuilt;
      for (const int n : in_cell) {
        particles[n].id() = -particles[n].id();
        ++n_invalid;
      }
      for (int k = 0; k < setup::n_lem; ++k) {
        ParticleType p;
        p.id()  = ParticleType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();
        for (int comp = 0; comp < RealData::ncomps; ++comp) {
          p.rdata(comp) = 0.0;
        }
        PlaceOnLine(p, iv, k, setup::n_lem, prob_lo, cell_size);
        p.idata(IntData::entry_side) = 0;
        fresh.push_back(p);
      }
    }
    // after the loop: push_back may move the storage of the tile
    for (const auto& p : fresh) {
      tile.push_back(p);
    }
  }

  amrex::ParallelDescriptor::ReduceLongSum(n_invalid);
  if (n_invalid > 0) {
    container.Redistribute(lev); // removes the replaced elements
  }
  container.RebuildCellMap(lev);
}

// element k of every cell of the fine level, as the sums above
void
PackElement(const ClemContainer& container, const int fine, const int k, amrex::MultiFab& packed)
{
  const auto& pc = container.Particles();

  for (ClemParConstIter pti(pc, fine); pti.isValid(); ++pti) {
    const amrex::Box& bx  = pti.validbox();
    const auto* particles = pti.GetArrayOfStructs().data();
    const auto& out       = packed.array(pti.index());

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      const auto& in_cell = container.ParticlesInCell(fine, pti.index(), iv);
      const int n         = FindElement(in_cell, particles, k);
      if (static_cast<int>(in_cell.size()) != setup::n_lem || n < 0) {
        amrex::AllPrint() << "CLEM average down: level " << fine << " cell " << iv << " holds "
                          << in_cell.size() << " elements and no regular line\n";
        amrex::Abort("CLEM average down: the fine level is not regridded");
      }
      const auto& p          = particles[n];
      const amrex::Real mass = p.rdata(RealData::mass);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        mass > 0.0 && std::abs(mass - particles[in_cell[0]].rdata(RealData::mass)) <=
                        equal_mass_tol * mass,
        "CLEM average down: a fine line is not of equal mass; regrid first");

      out(iv, c_mass)     = mass;
      out(iv, c_vol)      = p.rdata(RealData::vol);
      out(iv, c_mass_eint) = mass * p.rdata(RealData::eint);
      out(iv, c_mass_T)    = mass * p.rdata(RealData::T);
      out(iv, c_mass_phi)  = mass * p.rdata(RealData::phi);
      for (int s = 0; s < NUM_SPECIES; ++s) {
        out(iv, c_mass_Y + s) = mass * p.rdata(RealData::Y0 + s);
      }
    }
  }
}

// element k of every covered cell of level lev from the sums over its children
void
WriteElement(ClemContainer& container, const int lev, const int k, const amrex::MultiFab& sums,
             const amrex::iMultiFab& covered)
{
  auto& pc = container.Particles();
  auto eos = pele::physics::PhysicsType::eos();

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const amrex::Box& bx = pti.validbox();
    auto* particles      = pti.GetArrayOfStructs().data();
    const auto& sum      = sums.const_array(pti.index());
    const auto& mask     = covered.const_array(pti.index());

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (mask(iv) != 1) {
        continue;
      }
      const int n = FindElement(container.ParticlesInCell(lev, pti.index(), iv), particles, k);
      const amrex::Real mass = sum(iv, c_mass);
      const amrex::Real vol  = sum(iv, c_vol);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n >= 0 && mass > 0.0 && vol > 0.0,
                                       "CLEM average down: no fine line under a covered cell");

      auto& p                 = particles[n];
      p.rdata(RealData::mass) = mass;
      p.rdata(RealData::vol)  = vol;
      p.rdata(RealData::rho)  = mass / vol;
      p.rdata(RealData::eint) = sum(iv, c_mass_eint) / mass;
      p.rdata(RealData::phi)  = sum(iv, c_mass_phi) / mass;
      amrex::Real* Y          = &p.rdata(RealData::Y0);
      for (int s = 0; s < NUM_SPECIES; ++s) {
        Y[s] = sum(iv, c_mass_Y + s) / mass;
      }
      // the energy is the averaged one; the temperature follows from it
      p.rdata(RealData::T) = sum(iv, c_mass_T) / mass;
      TemperatureFromEnergy(p.rdata(RealData::rho), p.rdata(RealData::eint), Y,
                            p.rdata(RealData::T));
      eos.RTY2P(p.rdata(RealData::rho), p.rdata(RealData::T), Y, p.rdata(RealData::press));
    }
  }
}

} // namespace

AverageDiagnostics
AverageDown(ClemContainer& container, const int lev)
{
  BL_PROFILE("clem::AverageDown()");
  AMREX_ALWAYS_ASSERT(container.IsDefined());
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev >= 0 && lev < container.FinestLevel(),
                                   "AverageDown: lev must have a finer level above it");

  const int fine                 = lev + 1;
  const auto& pc                 = container.Particles();
  const amrex::BoxArray& crse_ba = pc.ParticleBoxArray(lev);
  const auto& crse_dm            = pc.ParticleDistributionMap(lev);
  const amrex::IntVect ref_ratio = pc.GetParGDB()->refRatio(lev);
  const amrex::Real n_children =
    static_cast<amrex::Real>(AMREX_D_TERM(ref_ratio[0], *ref_ratio[1], *ref_ratio[2]));

  // 1 in the cells of lev covered by the fine level
  const amrex::iMultiFab covered =
    amrex::makeFineMask(crse_ba, crse_dm, pc.ParticleBoxArray(fine), ref_ratio, 0, 1);

  AverageDiagnostics diag;
  PrepareCoveredCells(container, lev, covered, diag);
  if (!container.IsCellMapValid(fine)) {
    container.RebuildCellMap(fine);
  }

  // one element of the lines at a time, so the fields hold n_sums components
  // and not n_lem times as many
  amrex::MultiFab packed(pc.ParticleBoxArray(fine), pc.ParticleDistributionMap(fine), n_sums, 0);
  amrex::MultiFab sums(crse_ba, crse_dm, n_sums, 0);
  for (int k = 0; k < setup::n_lem; ++k) {
    PackElement(container, fine, k, packed);
    // the children on the rank of their parent: the mean times r^d is the sum
    sums.setVal(0.0);
    amrex::average_down(packed, sums, 0, n_sums, ref_ratio);
    sums.mult(n_children);
    WriteElement(container, lev, k, sums, covered);
  }

  amrex::ParallelDescriptor::ReduceLongSum(diag.n_covered);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_rebuilt);
  return diag;
}

void
AverageDownAll(ClemContainer& container)
{
  for (int lev = container.FinestLevel() - 1; lev >= 0; --lev) {
    AverageDown(container, lev);
  }
}

} // namespace clem
