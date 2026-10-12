/**
 * react.cpp: the reactor of the LEM lines.
 *
 *   ClemContainer::DiffusionStirringReactor  every line of a level advanced
 *                                            by dt
 *
 * The three operators that act ALONG a line are advanced together. Each line
 * is gathered into the workspace of the rank (LineScratch) and advanced by a
 * loop of explicit sub-steps that completes the LES dt:
 *
 *   while ( t_line < dt )
 *     dt_sub = min( dt_diffusion, dt - t_line )      [dt_stirring: to come]
 *     diffusion over dt_sub
 *     [stirring, reaction: to come]
 *     t_line += dt_sub
 *
 * and then written back to its particles. The models are those of setup.H.
 */
#include <AMReX_ParallelDescriptor.H>

#include "particle_container.H"
#include "line_scratch.H"
#include "setup.H"

namespace clem {

LineDiagnostics
ClemContainer::DiffusionStirringReactor(const int lev, const amrex::Real dt,
                                        const TransParmType* trans_parm)
{
  BL_PROFILE("clem::ClemContainer::DiffusionStirringReactor()");
  AMREX_ALWAYS_ASSERT(m_pc != nullptr);

  if (!IsCellMapValid(lev)) {
    RebuildCellMap(lev);
  }

  // the line of a cell runs along x: its length is the cell size, its
  // cross-section A = V_cell / length
  const auto cell_size       = m_pc->Geom(lev).CellSizeArray();
  const amrex::Real cell_vol = AMREX_D_TERM(cell_size[0], *cell_size[1], *cell_size[2]);
  const amrex::Real area     = cell_vol / cell_size[0];
  const amrex::Real area2    = area * area;

  const setup::DiffusionType diffusion;

  LineDiagnostics diag;
  auto& line = LineScratch::Get();

  amrex::Real enthalpy_before = 0.0;
  amrex::Real enthalpy_after  = 0.0;
  amrex::Real enthalpy_abs    = 0.0;
  amrex::Real unused          = 0.0;

  // the loop ends when less than this fraction of dt is left
  constexpr amrex::Real time_tol = 1.0e-12;

  for (ClemParIter pti(*m_pc, lev); pti.isValid(); ++pti) {
    const int grid_id    = pti.index();
    const amrex::Box& bx = pti.validbox();
    auto* particles      = pti.GetArrayOfStructs().data();

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (!line.Gather(ParticlesInCell(lev, grid_id, iv), particles)) {
        ++diag.n_skipped;
        continue;
      }
      line.Enthalpy(enthalpy_before, unused);

      amrex::Real t_line = 0.0;
      int n_sub          = 0;
      bool capped        = false;

      while (dt - t_line > time_tol * dt) {
        diffusion.UpdateProperties(line, trans_parm);

        // the sub-step: the smallest of the limits of the operators, and
        // never beyond the end of the LES step
        amrex::Real dt_sub          = diffusion.StableSubStep(line, area2);
        const amrex::Real remaining = dt - t_line;
        if (dt_sub >= remaining) {
          dt_sub = remaining;
        } else if (n_sub + 1 >= setup::max_substeps) {
          // out of sub-steps: finish the step above the stability limit
          dt_sub = remaining;
          capped = true;
        }

        diffusion.ComputeRates(line, area2);
        diffusion.AdvanceLine(line, dt_sub, diag);

        t_line += dt_sub;
        ++n_sub;
      }

      if (capped) {
        ++diag.n_capped;
      }
      diag.max_substeps = amrex::max<int>(diag.max_substeps, n_sub);
      diag.total_substeps += n_sub;
      ++diag.n_lines;

      line.Enthalpy(enthalpy_after, enthalpy_abs);
      line.Scatter(particles);
    }
  }

  amrex::ParallelDescriptor::ReduceLongSum(diag.n_lines);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_skipped);
  amrex::ParallelDescriptor::ReduceLongSum(diag.total_substeps);
  amrex::ParallelDescriptor::ReduceIntMax(diag.max_substeps);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_capped);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_clipped);
  amrex::ParallelDescriptor::ReduceRealMax(diag.max_neg_Y);
  amrex::ParallelDescriptor::ReduceRealSum(enthalpy_before);
  amrex::ParallelDescriptor::ReduceRealSum(enthalpy_after);
  amrex::ParallelDescriptor::ReduceRealSum(enthalpy_abs);
  diag.drift_enthalpy = (enthalpy_abs > 0.0) ? std::abs(enthalpy_after - enthalpy_before) / enthalpy_abs : 0.0;

  if (diag.n_capped > 0) {
    amrex::Warning("CLEM: lines hit setup::max_substeps and were advanced above their stability "
                   "limit; raise it or coarsen the line");
  }
  return diag;
}

} // namespace clem
