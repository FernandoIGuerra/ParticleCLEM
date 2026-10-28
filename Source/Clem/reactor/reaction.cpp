/**
 * reaction.cpp: the chemistry of the lines, through the CVODE of PelePhysics.
 *
 *   ClemContainer::InitializeReactor        the reactor of this rank, once
 *   ClemContainer::SetReactorTypicalValues  the typical values of its
 *                                           tolerances (PeleC's)
 *   DoReaction                              every element of a line over
 *                                           one sub-step
 *
 * The reactor is built as PeleC builds its own (ReactorBase::create, held by
 * pointer), but it is a second instance: the LES species are frozen, the
 * chemistry happens on the lines only. ONE per rank (ClemPC is one per rank),
 * and its CVODE context is built ONCE (ReactorCvode::InitClemReactor) and
 * reused by every element of every line of every sub-step through
 * CVodeReInit (ReactorCvode::ReactClem).
 */
#include "ReactorCvode.H"

#include "line_operators.H"

namespace clem {

using pele::physics::reactions::ReactorCvode;
using pele::physics::reactions::ReactorTypes;

// the element is advanced at constant volume: rho and e are held, Y_k and T
// advance, so nothing but the line changes after the reaction
static_assert(setup::reactor_type == ReactorTypes::e_reactor_type,
              "CLEM: only the constant-volume reactor (e_reactor_type) is supported");

ClemContainer::~ClemContainer()
{
  if (m_reactor != nullptr) {
    m_reactor->close();
  }
}

void
ClemContainer::InitializeReactor()
{
  BL_PROFILE("clem::ClemContainer::InitializeReactor()");
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_reactor == nullptr, "InitializeReactor: the reactor exists already");

  m_reactor = pele::physics::reactions::ReactorBase::create("ReactorCvode");
  m_cvode   = dynamic_cast<ReactorCvode*>(m_reactor.get());
  AMREX_ALWAYS_ASSERT(m_cvode != nullptr);

  // the ode.* and cvode.* inputs, the same as the reactor of PeleC
  m_reactor->init(setup::reactor_type, 1);
  m_cvode->InitClemReactor();
}

void
ClemContainer::SetReactorTypicalValues(const std::vector<amrex::Real>& typical_values)
{
  AMREX_ALWAYS_ASSERT(m_cvode != nullptr);
  m_reactor->set_typ_vals_ode(typical_values);
  m_cvode->SetClemTolerances();
}

void
DoReaction(ReactorCvode& reactor, LineScratch& line, const amrex::Real dt_sub, const amrex::Real time,
           LineDiagnostics& diag)
{
  for (int l = 0; l < line.n_elem; ++l) {
    // rho and eint held; Y and T of the element advanced in place
    const int rhs_evals = reactor.ReactClem<setup::reactor_type>(
      line.rho[l], &line.Y[l * NUM_SPECIES], line.T[l], line.eint[l], dt_sub, time);

    diag.rhs_evals += rhs_evals;
    diag.max_rhs_evals = amrex::max<int>(diag.max_rhs_evals, rhs_evals);
  }
  diag.n_reacted += line.n_elem;
}

} // namespace clem
