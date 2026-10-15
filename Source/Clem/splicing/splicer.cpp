#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_iMultiFab.H>

#include "PelePhysics.H"

#include "splicer.H"
#include "setup.H"

namespace clem {

namespace {

using ParticleType = ClemContainer::ParticleType;
// the line of one cell, head first
using Line = std::vector<ParticleType>;

// ---- faces of a cell: f = 2 d + side, side 0 the low face, 1 the high face ----
constexpr int n_face = 2 * AMREX_SPACEDIM;

int FaceDir(const int f) { return f / 2; }
// sign of the outward normal along +d
int FaceSign(const int f) { return (f % 2 == 0) ? -1 : 1; }
// the same face seen from the neighbour
int OppositeFace(const int f) { return f ^ 1; }

// entry_side of an element: 0 a resident, 1 + f arrived through face f of its cell
constexpr int resident = 0;
int EntryCode(const int f) { return 1 + f; }

amrex::IntVect
NeighbourCell(const amrex::IntVect& iv, const int f)
{
  amrex::IntVect iv_n = iv;
  iv_n[FaceDir(f)] += FaceSign(f);
  return iv_n;
}

// what lies behind a face: a cell of the level, or the boundary of the lines
// (0 on the valid cells and on the ghost cells covered by the level, 1 on the
// ghost cells that are not: outside the domain or on a coarser level)
constexpr int on_level  = 0;
constexpr int off_level = 1;

// outward transfer mass of each face of iv over the step: > 0 an outflux,
// < 0 an influx. flux[d] is signed along +d, on the face of index iv + side e_d
void
OutwardTransfers(const amrex::IntVect& iv, const amrex::Real dt,
                 const amrex::GpuArray<amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>& flux,
                 amrex::Real* mu_out)
{
  for (int f = 0; f < n_face; ++f) {
    const int d       = FaceDir(f);
    amrex::IntVect fv = iv;
    fv[d] += (f % 2); // [x_l x_r, y_l, y_r, z_l, z_r] = [0, 1, 2, 3, 4, 5]
    mu_out[f] = static_cast<amrex::Real>(FaceSign(f)) * dt * flux[d](fv);
  }
}

// the elements of one cell in lem_index order
void
GatherLine( const std::vector<int>& in_cell, 
            const ParticleType* particles, 
            Line& line)
{
  line.clear();
  for (const int n : in_cell) {
    line.push_back(particles[n]);
  }
  std::stable_sort(line.begin(), line.end(), [](const ParticleType& a, const ParticleType& b) {
    return a.idata(IntData::lem_index) < b.idata(IntData::lem_index);
  });
}

amrex::Real
LineMass(const Line& line)
{
  amrex::Real mass = 0.0;
  for (const auto& p : line) {
    mass += p.rdata(RealData::mass);
  }
  return mass;
}

// a piece of an element: the same state, its own mass and volume
ParticleType
Piece(const ParticleType& p, const amrex::Real mass)
{
  ParticleType piece         = p;
  piece.rdata(RealData::mass) = mass;
  piece.rdata(RealData::vol)  = mass / p.rdata(RealData::rho);
  return piece;
}

// Cut_mu(P) = (sigma, P^r): sigma the head of mass mu, P^r what is left
void
Cut(const Line& line, const amrex::Real mu, Line& sigma, Line& rest)
{
  sigma.clear();
  rest.clear();
  amrex::Real taken = 0.0;
  for (const auto& p : line) {
    const amrex::Real mass = p.rdata(RealData::mass);
    const amrex::Real want = mu - taken;
    if (want <= 0.0) {
      rest.push_back(p);
    } else if (mass <= want) {
      sigma.push_back(p);
      taken += mass;
    } else {
      sigma.push_back(Piece(p, want));
      rest.push_back(Piece(p, mass - want));
      taken = mu;
    }
  }
}

// the elements of a segment become particles of cell iv, in flight: at the
// centre of the cell, lem_index their order in the segment
void
Deliver(const Line& segment, const amrex::IntVect& iv, const int entry,
        std::vector<ParticleType>& out,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cell_size)
{
  for (int s = 0; s < static_cast<int>(segment.size()); ++s) {
    ParticleType p = segment[s];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      p.pos(d) = prob_lo[d] + (static_cast<amrex::Real>(iv[d]) + 0.5) * cell_size[d];
    }
    p.idata(IntData::lem_index)  = s;
    p.idata(IntData::entry_side) = entry;
    out.push_back(p);
  }
}

// one element of mass `mass` with the state of the LES cell iv_les
ParticleType
ElementFromLes(const amrex::Real mass, const amrex::Array4<const amrex::Real>& state,
               const amrex::IntVect& iv_les, const LesComps& comps)
{
  auto eos = pele::physics::PhysicsType::eos();

  const amrex::Real rho  = state(iv_les, comps.rho);
  const amrex::Real temp = state(iv_les, comps.temp);
  if (!(rho > 0.0 && temp > 0.0)) {
    amrex::Real rho_Y_sum = 0.0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      rho_Y_sum += state(iv_les, comps.spec + k);
    }
    amrex::AllPrint() << "CLEM splice: LES ghost cell " << iv_les << " gives " << mass
                      << " with rho " << rho << ", T " << temp << ", sum rho Y " << rho_Y_sum
                      << ", rho e " << state(iv_les, comps.eint) << '\n';
    amrex::Abort("CLEM splice: inflow from a LES ghost cell with no state");
  }

  ParticleType p;
  p.id()  = ParticleType::NextID();
  p.cpu() = amrex::ParallelDescriptor::MyProc();
  amrex::Real* Y = &p.rdata(RealData::Y0);
  for (int k = 0; k < NUM_SPECIES; ++k) {
    Y[k] = state(iv_les, comps.spec + k) / rho;
  }
  p.rdata(RealData::rho)  = rho;
  p.rdata(RealData::mass) = mass;
  p.rdata(RealData::vol)  = mass / rho;
  p.rdata(RealData::T)    = temp;
  p.rdata(RealData::phi)  = 0.0;
  eos.RTY2E(rho, temp, Y, p.rdata(RealData::eint));
  eos.RTY2P(rho, temp, Y, p.rdata(RealData::press));
  return p;
}

// the tile becomes exactly `out`
template <typename Tile>
void
ReplaceTile(Tile& tile, const std::vector<ParticleType>& out)
{
  tile.resize(out.size());
  auto& aos = tile.GetArrayOfStructs();
  for (std::size_t n = 0; n < out.size(); ++n) {
    aos[n] = out[n];
  }
}

// A^omega on every cell of the level: the tiles are rebuilt with what is left
// of each line, the segments addressed to their receptors and the inflow
// blocks. Local (unreduced) diagnostics
void
AdvectLines(ClemContainer& container, const int lev, const amrex::Real dt,
            const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& mass_flux,
            const amrex::MultiFab& les_border, const LesComps& comps, const int step,
            SpliceDiagnostics& diag)
{
  auto& pc                    = container.Particles();
  const amrex::Geometry& geom = pc.Geom(lev);
  const auto prob_lo          = geom.ProbLoArray();
  const auto cell_size        = geom.CellSizeArray();
  const setup::ReindexerType reindexer;

  // what lies behind each face, read in the ghost cells
  amrex::iMultiFab behind(pc.ParticleBoxArray(lev), pc.ParticleDistributionMap(lev), 1, 1);
  behind.BuildMask(geom.Domain(), geom.periodicity(), on_level, off_level, off_level, on_level);

  Line line;
  Line sigma;
  Line rest;
  Line next;
  std::vector<ParticleType> out;
  amrex::Real mu_out[n_face];
  int schedule[n_face];

  for (amrex::MFIter mfi = pc.MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const amrex::Box& bx  = mfi.validbox();
    const int grid_id     = mfi.index();
    auto& tile            = pc.DefineAndReturnParticleTile(lev, grid_id, mfi.LocalTileIndex());
    const auto* particles = tile.GetArrayOfStructs().data();
    const auto& state     = les_border.const_array(mfi);
    const auto& mask      = behind.const_array(mfi);
    
    const amrex::GpuArray<amrex::Array4<const amrex::Real>, AMREX_SPACEDIM> flux{AMREX_D_DECL(
      mass_flux[0].const_array(mfi), mass_flux[1].const_array(mfi), mass_flux[2].const_array(mfi))};

    out.clear();
    out.reserve(tile.numParticles() + n_face * bx.numPts());

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      //in line, the particle are collected in lem_index order, so that the cut is always at the head of the line
      GatherLine(container.ParticlesInCell(lev, grid_id, iv), particles, line);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!line.empty(), "CLEM splice: a cell has no line");
      const amrex::Real line_mass = LineMass(line);
      diag.mass_before            += line_mass;

      //computes from the face stored fluxes, the fluxx in each face \Gamma in a cell
      OutwardTransfers(iv, dt, flux, mu_out);

      // N^+_c in the order O^omega, and the donor cannot export more than it holds
      int n_out            = 0;
      amrex::Real mass_out = 0.0;
      for (int f = 0; f < n_face; ++f) {
        if (mu_out[f] > 0.0) {
          schedule[n_out++] = f;
          mass_out += mu_out[f];
        }
      }

      if (mass_out > line_mass * (1.0 + 1.0e-12)) {
        amrex::AllPrint() << "CLEM splice: level " << lev << " cell " << iv << " exports " << mass_out << " but holds " << line_mass << '\n';
        amrex::Abort("CLEM splice: a cell exports more than it holds; lower pelec.cfl");
      }

      diag.max_export = amrex::max<amrex::Real>(diag.max_export, mass_out / line_mass);
      reindexer.SortFlux(schedule, mu_out, n_out, step);

      // A^omega: Cut along the schedule, each segment to its neighbour
      rest = line;
      for (int l = 0; l < n_out; ++l) {
        const int f = schedule[l];
        Cut(rest, mu_out[f], sigma, next);
        rest.swap(next);

        const amrex::IntVect iv_n = NeighbourCell(iv, f);
        if (mask(iv_n) == on_level) {
          // a periodic image lies outside the domain: Redistribute wraps it
          Deliver(sigma, iv_n, EntryCode(OppositeFace(f)), out, prob_lo, cell_size);
        } else {
          diag.mass_out += LineMass(sigma); // the segment leaves the lines
        }
      }
      Deliver(rest, iv, resident, out, prob_lo, cell_size);

      // an influx with no line behind the face: one block from the LES ghost cell
      for (int f = 0; f < n_face; ++f) {
        const amrex::IntVect iv_n = NeighbourCell(iv, f);
        if (mu_out[f] >= 0.0 || mask(iv_n) == on_level) {
          continue;
        }
        line.assign(1, ElementFromLes(-mu_out[f], state, iv_n, comps));
        Deliver(line, iv, EntryCode(f), out, prob_lo, cell_size);
        diag.mass_in += -mu_out[f];
      }
    }
    ReplaceTile(tile, out);
  }
}

// R_I^omega on every cell, once the segments are delivered: what is left of
// the line, then the blocks ranked by O^omega(N^-_c); relabelled, and placed
// along x by their mass coordinate
void
ReindexLines(ClemContainer& container, const int lev, const amrex::Real dt,
             const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& mass_flux, const int step,
             SpliceDiagnostics& diag)
{
  auto& pc             = container.Particles();
  const auto prob_lo   = pc.Geom(lev).ProbLoArray();
  const auto cell_size = pc.Geom(lev).CellSizeArray();
  const setup::ReindexerType reindexer;

  std::vector<int> line;
  amrex::Real mu_out[n_face];
  amrex::Real mu_in[n_face];
  int order[n_face];
  int rank[n_face];

  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const amrex::Box& bx = pti.validbox();
    auto* particles      = pti.GetArrayOfStructs().data();
    const amrex::GpuArray<amrex::Array4<const amrex::Real>, AMREX_SPACEDIM> flux{AMREX_D_DECL(
      mass_flux[0].const_array(pti), mass_flux[1].const_array(pti), mass_flux[2].const_array(pti))};

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      line = container.ParticlesInCell(lev, pti.index(), iv);

      // O^omega on N^-_c: rank 1, 2, ... of each influx face; 0 is the residual
      OutwardTransfers(iv, dt, flux, mu_out);
      int n_in = 0;
      for (int f = 0; f < n_face; ++f) {
        mu_in[f] = amrex::max<amrex::Real>(-mu_out[f], 0.0);
        rank[f]  = 0;
        if (mu_in[f] > 0.0) {
          order[n_in++] = f;
        }
      }
      reindexer.SortIncomingParticle(order, mu_in, n_in, step);
      for (int l = 0; l < n_in; ++l) {
        rank[order[l]] = l + 1;
      }

      auto rank_of = [&](const int n) {
        const int entry = particles[n].idata(IntData::entry_side);
        return (entry == resident) ? 0 : rank[entry - 1];
      };
      
      amrex::Real line_mass = 0.0;
      for (const int n : line) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          particles[n].idata(IntData::entry_side) == resident || rank_of(n) > 0,
          "CLEM splice: a block arrived through a face with no influx");
        line_mass += particles[n].rdata(RealData::mass);
      }
      std::stable_sort(line.begin(), line.end(), [&](const int a, const int b) {
        const int rank_a = rank_of(a);
        const int rank_b = rank_of(b);
        if (rank_a != rank_b) {
          return rank_a < rank_b;
        }
        return particles[a].idata(IntData::lem_index) < particles[b].idata(IntData::lem_index);
      });

      // lem_index follows from the position; nothing is in flight any more
      amrex::Real mass_before = 0.0;
      for (int q = 0; q < static_cast<int>(line.size()); ++q) {
        auto& p                      = particles[line[q]];
        p.idata(IntData::lem_index)  = q;
        p.idata(IntData::entry_side) = resident;
        // strictly inside the cell, so AMReX keeps it there
        const amrex::Real frac =
          (line_mass > 0.0) ? (mass_before + 0.5 * p.rdata(RealData::mass)) / line_mass : 0.5;
        p.pos(0) = prob_lo[0] + (static_cast<amrex::Real>(iv[0]) + 0.05 + 0.9 * frac) * cell_size[0];
        mass_before += p.rdata(RealData::mass);
      }
      diag.mass_after += line_mass;
      diag.n_elements += static_cast<amrex::Long>(line.size());
    }
  }
}

} // namespace

SpliceDiagnostics
Splice(ClemContainer& container, const int lev, const amrex::Real dt,
       const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& mass_flux,
       const amrex::MultiFab& les_border, const LesComps& comps, const int step)
{
  BL_PROFILE("clem::Splice()");
  const auto& pc = container.Particles();
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    les_border.boxArray() == pc.ParticleBoxArray(lev) &&
      les_border.DistributionMap() == pc.ParticleDistributionMap(lev) && les_border.nGrow() >= 1,
    "Splice: LES state must live on the particle grids of this level, with a ghost cell");
  
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      mass_flux[d].boxArray() == amrex::convert(pc.ParticleBoxArray(lev), amrex::IntVect::TheDimensionVector(d)), 
      "Splice: mass_flux must live on the faces of the particle grids of this level");
  }

  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev); // which element lies in which cell
  }
  SpliceDiagnostics diag;
  AdvectLines(container, lev, dt, mass_flux, les_border, comps, step, diag); // A^omega
  container.Redistribute(lev);   // delivers every segment to its receptor
  container.RebuildCellMap(lev); // the lines with their arrivals

  ReindexLines(container, lev, dt, mass_flux, step, diag); // R_I^omega

  amrex::Real sums[4] = {diag.mass_before, diag.mass_after, diag.mass_in, diag.mass_out};
  amrex::ParallelDescriptor::ReduceRealSum(sums, 4);
  diag.mass_before = sums[0];
  diag.mass_after  = sums[1];
  diag.mass_in     = sums[2];
  diag.mass_out    = sums[3];
  amrex::ParallelDescriptor::ReduceRealMax(diag.max_export);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_elements);
  return diag;
}

amrex::Real
LineMassDefect(ClemContainer& container, const int lev, const amrex::MultiFab& les,
               const int rho_comp)
{
  auto& pc = container.Particles();
  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev);
  }
  const auto cell_size       = pc.Geom(lev).CellSizeArray();
  const amrex::Real cell_vol = AMREX_D_TERM(cell_size[0], *cell_size[1], *cell_size[2]);

  amrex::Real max_defect = 0.0;
  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const amrex::Box& bx  = pti.validbox();
    const auto* particles = pti.GetArrayOfStructs().data();
    const auto& state     = les.const_array(pti);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      amrex::Real line_mass = 0.0;
      for (const int n : container.ParticlesInCell(lev, pti.index(), iv)) {
        line_mass += particles[n].rdata(RealData::mass);
      }
      const amrex::Real cell_mass = state(iv, rho_comp) * cell_vol;
      max_defect = amrex::max<amrex::Real>(max_defect, std::abs(line_mass / cell_mass - 1.0));
    }
  }
  amrex::ParallelDescriptor::ReduceRealMax(max_defect);
  return max_defect;
}

amrex::Real
LineTemperatureGap(ClemContainer& container, const int lev, const amrex::MultiFab& les,
                   const int temp_comp)
{
  auto& pc = container.Particles();
  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev);
  }

  amrex::Real max_gap = 0.0;
  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const amrex::Box& bx  = pti.validbox();
    const auto* particles = pti.GetArrayOfStructs().data();
    const auto& state     = les.const_array(pti);

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      amrex::Real mass   = 0.0;
      amrex::Real mass_T = 0.0;
      for (const int n : container.ParticlesInCell(lev, pti.index(), iv)) {
        mass += particles[n].rdata(RealData::mass);
        mass_T += particles[n].rdata(RealData::mass) * particles[n].rdata(RealData::T);
      }
      if (mass > 0.0) {
        max_gap = amrex::max<amrex::Real>(max_gap, std::abs(mass_T / mass - state(iv, temp_comp)));
      }
    }
  }
  amrex::ParallelDescriptor::ReduceRealMax(max_gap);
  return max_gap;
}

} // namespace clem
