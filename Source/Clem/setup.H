/**
 * Compile-time setup of the CLEM model: every choice that defines what a
 * build evaluates lives here.
 */
#ifndef CLEM_SETUP_H
#define CLEM_SETUP_H

#include <AMReX_REAL.H>

#include "ReactorTypes.H"

#include "mixture_averaged_diffusion.H"
#include "write_back_coupling.H"
#include "equal_mass_regridder.H"
#include "positive_reindexer.H"
#include "negative_reindexer.H"
#include "alternating_reindexer.H"

namespace clem::setup {

// number of LEM elements on the line of one LES cell
static constexpr int n_lem = 12;

// ---- coupling between the LES and the lines (coupling/) -----------------------

// WriteBack: isentropic pressure forward, conserved-delta write-back
enum class Coupling { WriteBack };

static constexpr Coupling coupling = Coupling::WriteBack;

template <Coupling C>
struct CouplingOf;
template <>
struct CouplingOf<Coupling::WriteBack>
{
  using type = clem::WriteBackCoupling;
};

using CouplingType = typename CouplingOf<coupling>::type;

// return the lines to the LES at the end of the CLEM step (BackwardCoupling):
// rho Y_k and rho e of the LES become what the lines transported and
// processed. With ForwardCoupling this is the two-way coupling of the model
static constexpr bool couple_back = true;

// the LES does not advance rho Y_k: the species part of the MOL source is
// zeroed in both RK2 stages, so the species are transported by the lines
// alone and reach the LES only through BackwardCoupling. Needs couple_back
static constexpr bool freeze_les_species = true;
static_assert(!freeze_les_species || couple_back,
              "CLEM: freeze_les_species needs couple_back, or the LES species never change");

// the energy channel of BackwardCoupling; rho E is always rebuilt with the
// kinetic energy of the LES, whose momentum the lines do not touch
//   LesPlusSubgrid: (rho e)^{n+1} = (rho E** - k**) + Delta[rho e]_subgrid.
//                   The LES keeps what its energy equation did (transport,
//                   p dV work, viscous heating, conduction); the lines add only
//                   what the LES does not have: the delta is taken around the
//                   line processes alone (diffusion; stirring and reaction to
//                   come), not around the splice or the pressure coupling
//   LinesOnly:      (rho e)^{n+1} = (rho E^n - k^n) + Delta[rho e] of the whole
//                   CLEM step: the energy of the cell is that of the lines
//                   (CLEMNew, clem.couple_back_energy = 1). Unstable at a
//                   resolved jet lip, see multiresolution/EnergyCouplingProblem.pdf
//   Off:            rho E of the LES, no energy from the lines (diagnostic)
enum class EnergyCoupling { LesPlusSubgrid, LinesOnly, Off };

static constexpr EnergyCoupling energy_coupling = EnergyCoupling::LesPlusSubgrid;

// after the write-back, shift the internal energy of every element of a line
// by the same amount so that the line's mass-mean e is the e of its LES cell
// (AnchorLineEnergy). Closes the energy loop LES -> lines: without it the
// lines keep only isentropic pressure work and drift from the LES wherever
// the LES heats or cools irreversibly (shocks, walls, conduction)
static constexpr bool anchor_line_energy = true;

// ---- regridding of the lines (regridding/) ------------------------------------

// EqualMass: n_lem elements of equal mass; e, Y mass-averaged, T, p from the EOS
enum class Regridder { EqualMass };

static constexpr Regridder regridder = Regridder::EqualMass;

template <Regridder R>
struct RegridderOf;
template <>
struct RegridderOf<Regridder::EqualMass>
{
  using type = clem::EqualMassRegridder;
};

using RegridderType = typename RegridderOf<regridder>::type;

// ---- ordering policy of the splicing (reindexing/) -----------------------------

// Positive:    Sigma^+, the largest transfer first
// Negative:    Sigma^-, the smallest transfer first
// Alternating: Sigma^- on the even steps, Sigma^+ on the odd ones
enum class Reindexer { Positive, Negative, Alternating };

static constexpr Reindexer reindexer = Reindexer::Positive;

template <Reindexer R>
struct ReindexerOf;
template <>
struct ReindexerOf<Reindexer::Positive>
{
  using type = clem::PositiveReindexer;
};
template <>
struct ReindexerOf<Reindexer::Negative>
{
  using type = clem::NegativeReindexer;
};
template <>
struct ReindexerOf<Reindexer::Alternating>
{
  using type = clem::AlternatingReindexer;
};

using ReindexerType = typename ReindexerOf<reindexer>::type;

// ---- molecular diffusion along the line (diffusion/) ---------------------------

// MixtureAveraged: Y_k and T, mass coordinate, isobaric, periodic line
enum class Diffusion { MixtureAveraged };

static constexpr Diffusion diffusion = Diffusion::MixtureAveraged;

// the model type the rest of the code uses; a new model adds one enum value
// and one branch here
template <Diffusion D>
struct DiffusionOf;
template <>
struct DiffusionOf<Diffusion::MixtureAveraged>
{
  using type = clem::MixtureAveragedDiffusion;
};

using DiffusionType = typename DiffusionOf<diffusion>::type;

// diffusion number of the explicit sub-step, in (0, 0.5]: 0.5 is the
// stability limit of the scheme
static constexpr amrex::Real diffusion_cfl = 0.5;

// ---- operators along the lines (reactor/) ------------------------------------------

// largest number of sub-steps of one line in one LES step; a line that needs
// more finishes the step above its stability limit, and is reported
static constexpr int max_substeps = 1000;

// chemistry of the lines (reactor/reaction.cpp): every element is a 0-D
// reactor integrated by the CVODE of PelePhysics, over each sub-step
static constexpr bool react_lines = true;

// constant volume: rho and e of the element are held, Y_k and T advance
static constexpr int reactor_type = pele::physics::reactions::ReactorTypes::e_reactor_type;

/** @brief skip the reaction of the lines of covered cells, to reduce the
 *  cost: the average down replaces them at the sync. They still diffuse and
 *  splice */
static constexpr bool skip_covered = false;

// ---- LES closure (PeleC::read_params) -------------------------------------------

// The subgrid closure of the LES. In a CLEM build these choices REPLACE the LES
// inputs of the input file (pelec.do_les, les_model, Cw, CI, PrT,
// use_explicit_filter), whatever it says. PeleC has LES in 3D only: a 2D build
// runs without it, nu_T = 0 and the lines are not stirred
//   Wale: wall-adapting local eddy viscosity (Nicoud & Ducros 1999),
//         mu_T = rho (Cw Delta)^2 OP_WALE, Delta = (dx dy dz)^(1/3)
// the LES closure on or off; off, nu_T = 0 and the lines are not stirred
static constexpr bool do_les = true;

enum class LesModel { Wale };

static constexpr LesModel les_model = LesModel::Wale;

// the number of each model in PeleC (PeleC::getLESTerm)
constexpr int
PeleCLesModel(const LesModel model)
{
  return model == LesModel::Wale ? 2 : -1;
}

// WALE constant Cw
static constexpr amrex::Real wale_cw = 0.5;

// coefficient of the isotropic part of the subgrid stress, C_I; 0 leaves it out
static constexpr amrex::Real les_ci = 0.0;

// turbulent Prandtl number of the subgrid heat flux of the LES
static constexpr amrex::Real les_prt = 1.0;

// the LES is filtered implicitly, by its grid
static constexpr bool les_explicit_filter = false;

// turbulent Schmidt number of the lines: the turbulent diffusivity the
// stirring reproduces is D_T = nu_T / Sc_T, with nu_T = mu_T / rho of the
// LES closure (PeleC::clem_nu_t)
static constexpr amrex::Real turbulent_schmidt = 1.0;

// stirring of the lines: triplet maps at the rate that reproduces the D_T of
// the LES closure (stirring/, reactor/stirring.cpp)
static constexpr bool stir_lines = true;

// the fewest elements a triplet map acts on (a map of 3 is the identity):
// eddies smaller than this many elements of a line, l < 6 Delta_x / n_lem,
// are left out of the inertial range of the stirring (stirring/eddy_sampling.H)
static constexpr int stirring_min_elements = 6;

// ---- checks ---------------------------------------------------------------------

// after every LES step, print how well the stored mass flux closes the LES
// mass balance of each level (PeleC::CheckClemMassFlux), and how well the
// mass of the lines follows the mass of their cells after the splice
static constexpr bool check_mass_flux = true;

} // namespace clem::setup
#endif
