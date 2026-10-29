/**
 * The operators that act ALONG one gathered LEM line, over one sub-step of
 * ClemContainer::DiffusionStirringReactor (diffusion_reaction_stirring.cpp).
 *
 * Functions:
 *
 * DoDiffusion : diffusion of the line over dt_sub          (diffusion.cpp)
 *
 * DoReaction  : chemistry of every element over dt_sub      (reaction.cpp)
 *
 * DoStirring  : the triplet map of one eddy on the line     (stirring.cpp)
 */
#ifndef CLEM_LINE_OPERATORS_H
#define CLEM_LINE_OPERATORS_H

#include "line_scratch.H"
#include "particle_container.H"
#include "setup.H"

namespace clem {

// DoDiffusion: the rates of the line (ComputeRates) and the update of its
// elements (AdvanceLine) over dt_sub. Needs diffusion.UpdateProperties first,
// the same that sized dt_sub
void DoDiffusion(const setup::DiffusionType& diffusion, LineScratch& line, amrex::Real dt_sub,
                 amrex::Real area2, LineDiagnostics& diag);

// DoReaction: every element of the line is a 0-D reactor of setup::reactor_type,
// integrated alone by the CVODE of this rank over [time, time + dt_sub]
void DoReaction(pele::physics::reactions::ReactorCvode& reactor, LineScratch& line,
                amrex::Real dt_sub, amrex::Real time, LineDiagnostics& diag);

// DoStirring: the triplet map of an eddy of n_eddy elements (a multiple of 3,
// at least setup::stirring_min_elements, at most the line) whose first
// element is `start`, on the periodic line. The eddy is drawn by the caller
// (stirring/eddy_sampling.H)
void DoStirring(LineScratch& line, int start, int n_eddy);

} // namespace clem
#endif
