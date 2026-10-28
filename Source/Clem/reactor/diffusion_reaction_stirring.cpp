/**
 * diffusion_reaction_stirring.cpp: the operators along the LEM lines.
 *
 *   ClemContainer::DiffusionStirringReactor  every line of a level advanced
 *                                            by dt
 *
 * The three operators that act ALONG a line are advanced together. Each line
 * is gathered into the workspace of the rank (LineScratch) and advanced by a
 * loop of explicit sub-steps that completes the LES dt:
 *
 *   while ( t_line < dt )
 *     dt_sub = min( dt_diffusion, t_stirring - t_line, dt - t_line )
 *     DoDiffusion over dt_sub                               (diffusion.cpp)
 *     DoReaction  over the same dt_sub                       (reaction.cpp)
 *     t_line += dt_sub
 *     if ( t_line >= t_stirring )  DoStirring                (stirring.cpp)
 *
 * and then written back to its particles. The models are those of setup.H.
 * Stirring is not implemented yet: t_stirring is never reached.
 */
#include <limits>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include "line_operators.H"

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

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!setup::react_lines || HasReactor(),
                                   "DiffusionStirringReactor: the reactor is not initialized");

  // 1 where the cell is covered by level lev + 1: its line is not reacted
  // if setup::skip_covered. The finest level covers nothing
  const bool has_covered = setup::skip_covered && lev < m_pc->finestLevel();
  amrex::iMultiFab covered;
  if (has_covered) {
    covered = amrex::makeFineMask(m_pc->ParticleBoxArray(lev), m_pc->ParticleDistributionMap(lev),
                                  m_pc->ParticleBoxArray(lev + 1), m_pc->GetParGDB()->refRatio(lev), 0, 1);
  }

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
    const auto mask      = has_covered ? covered.const_array(grid_id) : amrex::Array4<const int>{};

    for (amrex::IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
      if (!line.Gather(ParticlesInCell(lev, grid_id, iv), particles)) {
        ++diag.n_skipped;
        continue;
      }
      line.Enthalpy(enthalpy_before, unused);

      const bool react_line = setup::react_lines && !(has_covered && mask(iv) == 1);
      if (setup::react_lines && !react_line) {
        ++diag.n_skip_react;
      }

      amrex::Real t_line = 0.0;
      int n_sub          = 0;
      bool capped        = false;

      // time of the next stirring event of the line: none until stirring exists
      const amrex::Real t_stirring = std::numeric_limits<amrex::Real>::max();

      while (dt - t_line > time_tol * dt) {
        diffusion.UpdateProperties(line, trans_parm);

        // the sub-step: the smallest of the limits of the operators, and
        // never beyond the end of the LES step
        amrex::Real dt_sub          = amrex::min(diffusion.StableSubStep(line, area2), t_stirring - t_line);
        const amrex::Real remaining = dt - t_line;
        if (dt_sub >= remaining) {
          dt_sub = remaining;
        } else if (n_sub + 1 >= setup::max_substeps) {
          // out of sub-steps: finish the step above the stability limit
          dt_sub = remaining;
          capped = true;
        }

        DoDiffusion(diffusion, line, dt_sub, area2, diag);
        if (react_line) {
          DoReaction(*m_cvode, line, dt_sub, t_line, diag);
        }

        t_line += dt_sub;
        ++n_sub;

        if (t_line >= t_stirring) {
          DoStirring(line);
        }
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
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_reacted);
  amrex::ParallelDescriptor::ReduceLongSum(diag.n_skip_react);
  amrex::ParallelDescriptor::ReduceLongSum(diag.rhs_evals);
  amrex::ParallelDescriptor::ReduceIntMax(diag.max_rhs_evals);
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
