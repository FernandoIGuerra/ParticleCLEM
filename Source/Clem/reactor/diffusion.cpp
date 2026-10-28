/**
 * diffusion.cpp: the diffusion of one line over one sub-step.
 *
 *   DoDiffusion  ComputeRates and AdvanceLine of the model of setup.H, in one
 *                call
 */
#include "line_operators.H"

namespace clem {

void
DoDiffusion(const setup::DiffusionType& diffusion, LineScratch& line, const amrex::Real dt_sub,
            const amrex::Real area2, LineDiagnostics& diag)
{
  diffusion.ComputeRates(line, area2);
  diffusion.AdvanceLine(line, dt_sub, diag);
}

} // namespace clem
