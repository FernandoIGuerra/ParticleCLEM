#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_ParallelDescriptor.H>

#include "refluxer.H"
#include "setup.H"

namespace clem {

namespace {

using ParticleType = ClemContainer::ParticleType;

// faces of a cell: f = 2 d + side, side 0 the low face, 1 the high face
constexpr int n_face = 2 * AMREX_SPACEDIM;

int FaceDir(const int f) { return f / 2; }
int FaceSign(const int f) { return (f % 2 == 0) ? -1 : 1; }

amrex::IntVect
NeighbourCell(const amrex::IntVect& iv, const int f)
{
  amrex::IntVect iv_n = iv;
  iv_n[FaceDir(f)] += FaceSign(f);
  return iv_n;
}

// the transfer still to be made through each interface face of iv, outward:
// > 0, iv must give that mass; < 0, iv must receive it; 0 on the other faces
bool
OutwardMismatch(const amrex::IntVect& iv, const amrex::Array4<const int>& covered,
                const amrex::GpuArray<amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>& mismatch,
                amrex::Real* outward)
{
  bool any = false;
  for (int f = 0; f < n_face; ++f) {
    outward[f] = 0.0;
    if (covered(NeighbourCell(iv, f)) != 1) {
      continue;
    }
    const int d       = FaceDir(f);
    amrex::IntVect fv = iv;
    fv[d] += (f % 2);
    outward[f] = static_cast<amrex::Real>(FaceSign(f)) * mismatch[d](fv);
    any        = any || (outward[f] != 0.0);
  }
  return any;
}

// rank 1, 2, ... of the faces through which iv receives, in the order of the
// policy; 0 for the others
void
ArrivalRanks(const amrex::Real* outward, const int step, int* rank)
{
  const setup::ReindexerType reindexer;
  amrex::Real mu_in[n_face];
  int order[n_face];
  int n_in = 0;
  for (int f = 0; f < n_face; ++f) {
    mu_in[f] = amrex::max<amrex::Real>(-outward[f], 0.0);
    rank[f]  = 0;
    if (mu_in[f] > 0.0) {
      order[n_in++] = f;
    }
  }
  reindexer.SortIncomingParticle(order, mu_in, n_in, step);
  for (int l = 0; l < n_in; ++l) {
    rank[order[l]] = l + 1;
  }
}

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

} // namespace

RefluxDiagnostics
RefluxMass(ClemContainer& container, const int lev,
           const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& mismatch,
           const amrex::iMultiFab& covered, const int step, amrex::iMultiFab& touched)
{
  BL_PROFILE("clem::RefluxMass()");
  AMREX_ALWAYS_ASSERT(container.IsDefined());
  AMREX_ALWAYS_ASSERT(covered.nGrow() >= 1);

  auto& pc                       = container.Particles();
  const amrex::Geometry& geom    = pc.Geom(lev);
  const amrex::BoxArray& ba      = pc.ParticleBoxArray(lev);
  const auto& dm                 = pc.ParticleDistributionMap(lev);
  const auto prob_lo             = geom.ProbLoArray();
  const auto cell_size           = geom.CellSizeArray();

  if (!container.IsCellMapValid(lev)) {
    container.RebuildCellMap(lev);
  }

  RefluxDiagnostics diag;
  touched.setVal(0);

  // what each cell still has to receive through each of its faces
  amrex::MultiFab need(ba, dm, n_face, 0);
  need.setVal(0.0);

  amrex::Real outward[n_face];
  int rank[n_face];
  std::vector<int> line;
  amrex::Long n_invalid = 0;

  // ---- C gave too little: cut from its head, and dropped ----------------------
  for (ClemParIter pti(pc, lev); pti.isValid(); ++pti) {
    const amrex::Box& bx = pti.validbox();
    auto* particles      = pti.GetArrayOfStructs().data();
    const auto& mask     = covered.const_array(pti.index());
    const auto& mark     = touched.array(pti.index());
    const auto& want     = need.array(pti.index());
    const amrex::GpuArray<amrex::Array4<const amrex::Real>, AMREX_SPACEDIM> delta{AMREX_D_DECL(
      mismatch[0].const_array(pti), mismatch[1].const_array(pti), mismatch[2].const_array(pti))};

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (mask(iv) == 1 || !OutwardMismatch(iv, mask, delta, outward)) {
        continue;
      }
      mark(iv) = 1;
      ++diag.n_cells;

      line = container.ParticlesInCell(lev, pti.index(), iv);
      std::sort(line.begin(), line.end(), [&](const int a, const int b) {
        return particles[a].idata(IntData::lem_index) < particles[b].idata(IntData::lem_index);
      });
      amrex::Real line_mass = 0.0;
      for (const int n : line) {
        line_mass += particles[n].rdata(RealData::mass);
      }

      amrex::Real cut  = 0.0;
      amrex::Real gain = 0.0;
      for (int f = 0; f < n_face; ++f) {
        cut += amrex::max<amrex::Real>(outward[f], 0.0);
        want(iv, f) = amrex::max<amrex::Real>(-outward[f], 0.0);
        gain += want(iv, f);
      }
      if (cut > line_mass) {
        amrex::AllPrint() << "CLEM reflux: level " << lev << " cell " << iv << " must give "
                          << cut << " but holds " << line_mass << '\n';
        amrex::Abort("CLEM reflux: a mass mismatch is larger than the line it comes from");
      }
      diag.mass_removed += cut;
      diag.max_fraction =
        amrex::max<amrex::Real>(diag.max_fraction, std::abs(gain - cut) / line_mass);

      // all the cuts come from the head and all are dropped, so their order
      // does not matter: one cut of their sum
      for (const int n : line) {
        if (cut <= 0.0) {
          break;
        }
        auto& p                = particles[n];
        const amrex::Real mass = p.rdata(RealData::mass);
        if (mass <= cut) {
          p.id() = -p.id();
          ++n_invalid;
          cut -= mass;
        } else {
          p.rdata(RealData::mass) = mass - cut;
          p.rdata(RealData::vol)  = (mass - cut) / p.rdata(RealData::rho);
          cut                     = 0.0;
        }
      }
    }
  }

  // ---- C needs more: from the head of the covered cell behind the face --------
  // one element of the covered lines at a time, from the head: element k of
  // every covered cell as a field with a ghost cell, so that C reads the line
  // of its neighbour across grids and ranks. Most mismatches fit in element 0
  amrex::MultiFab donor(ba, dm, RealData::ncomps, 1);
  std::vector<ParticleType> fresh;
  for (int k = 0; k < setup::n_lem; ++k) {
    amrex::Real left = need.norm0(0, 0, false, true);
    for (int f = 1; f < n_face; ++f) {
      left = amrex::max<amrex::Real>(left, need.norm0(f, 0, false, true));
    }
    amrex::ParallelDescriptor::ReduceRealMax(left);
    if (left <= 0.0) {
      break;
    }

    donor.setVal(0.0);
    for (ClemParConstIter pti(pc, lev); pti.isValid(); ++pti) {
      const amrex::Box& bx  = pti.validbox();
      const auto* particles = pti.GetArrayOfStructs().data();
      const auto& mask      = covered.const_array(pti.index());
      const auto& out       = donor.array(pti.index());
      for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
        if (mask(iv) != 1) {
          continue;
        }
        const int n = FindElement(container.ParticlesInCell(lev, pti.index(), iv), particles, k);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n >= 0,
                                         "CLEM reflux: a covered cell has no regular line");
        for (int comp = 0; comp < RealData::ncomps; ++comp) {
          out(iv, comp) = particles[n].rdata(comp);
        }
      }
    }
    donor.FillBoundary(geom.periodicity());

    for (amrex::MFIter mfi = pc.MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.validbox();
      auto& tile           = pc.DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
      const auto& mask     = covered.const_array(mfi);
      const auto& mark     = touched.const_array(mfi);
      const auto& want     = need.array(mfi);
      const auto& src      = donor.const_array(mfi);
      const amrex::GpuArray<amrex::Array4<const amrex::Real>, AMREX_SPACEDIM> delta{AMREX_D_DECL(
        mismatch[0].const_array(mfi), mismatch[1].const_array(mfi), mismatch[2].const_array(mfi))};

      fresh.clear();
      for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
        if (mark(iv) != 1) {
          continue;
        }
        OutwardMismatch(iv, mask, delta, outward);
        ArrivalRanks(outward, step, rank);
        for (int f = 0; f < n_face; ++f) {
          if (want(iv, f) <= 0.0) {
            continue;
          }
          const amrex::IntVect iv_n = NeighbourCell(iv, f);
          const amrex::Real take =
            amrex::min<amrex::Real>(want(iv, f), src(iv_n, RealData::mass));
          if (take <= 0.0) {
            continue;
          }
          ParticleType p;
          p.id()  = ParticleType::NextID();
          p.cpu() = amrex::ParallelDescriptor::MyProc();
          for (int comp = 0; comp < RealData::ncomps; ++comp) {
            p.rdata(comp) = src(iv_n, comp);
          }
          p.rdata(RealData::mass) = take;
          p.rdata(RealData::vol)  = take / p.rdata(RealData::rho);
          PlaceOnLine(p, iv, 0, setup::n_lem, prob_lo, cell_size);
          // at the tail: after the residents, block by block in the order of
          // the policy, element by element inside a block
          p.idata(IntData::lem_index)  = setup::n_lem * rank[f] + k;
          p.idata(IntData::entry_side) = 0;
          fresh.push_back(p);

          want(iv, f) -= take;
          diag.mass_added += take;
        }
      }
      // after the loop: push_back may move the storage of the tile
      for (const auto& p : fresh) {
        tile.push_back(p);
      }
    }
  }

  amrex::Real left = need.norm0(0, 0, false, true);
  for (int f = 1; f < n_face; ++f) {
    left = amrex::max<amrex::Real>(left, need.norm0(f, 0, false, true));
  }
  amrex::ParallelDescriptor::ReduceRealMax(left);
  if (left > 0.0) {
    amrex::Abort("CLEM reflux: a mass mismatch is larger than the line it comes from");
  }

  amrex::ParallelDescriptor::ReduceLongSum(n_invalid);
  if (n_invalid > 0) {
    container.Redistribute(lev); // removes the elements cut whole
  }
  container.RebuildCellMap(lev);

  amrex::ParallelDescriptor::ReduceLongSum(diag.n_cells);
  amrex::ParallelDescriptor::ReduceRealSum(diag.mass_removed);
  amrex::ParallelDescriptor::ReduceRealSum(diag.mass_added);
  amrex::ParallelDescriptor::ReduceRealMax(diag.max_fraction);
  return diag;
}

} // namespace clem
