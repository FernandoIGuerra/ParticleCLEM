/**
 * Here we create the base class of the regridding of the LEM lines.
 *
 * Regridder<Derived>: the base of every regridding strategy (CRTP, no virtual
 * calls). The operators of the CLEM step leave a line with any number of
 * elements, of any mass (the splicing cuts and appends, the pressure and the
 * diffusion change the volumes). The regridding puts every line of ONE level
 * back on the grid the line operators assume:
 *
 *   Regrid  every cell of the level back to setup::n_lem elements
 *
 * The strategy is chosen at compile time in setup.H.
 */
#ifndef CLEM_REGRIDDER_H
#define CLEM_REGRIDDER_H

#include <AMReX_iMultiFab.H>

#include "particle_container.H"

namespace clem {

// what one Regrid did, over the level and the ranks
struct RegridDiagnostics
{
  amrex::Long n_added   = 0; // empty elements created in the cells that held too few
  amrex::Long n_removed = 0; // elements removed from the cells that held too many
};

template <typename DerivedRegridder>
class Regridder
{
  friend DerivedRegridder;

  // only a derived strategy can build the base
  Regridder() = default;

  const DerivedRegridder& derived() const { return static_cast<const DerivedRegridder&>(*this); }

public:
  ~Regridder()                           = default;
  Regridder(Regridder const&)            = delete;
  Regridder& operator=(Regridder const&) = delete;
  Regridder(Regridder&&)                 = delete;
  Regridder& operator=(Regridder&&)      = delete;

  // on return every cell of level lev holds exactly setup::n_lem elements and
  // the cell map of lev is valid. If only_cells is given, only the cells
  // where it is 1 are regridded and the others are left as they are
  RegridDiagnostics Regrid(ClemContainer& container, const int lev,
                           const amrex::iMultiFab* only_cells = nullptr) const
  {
    return derived().RegridImplementation(container, lev, only_cells);
  }

  static constexpr const char* Name() { return DerivedRegridder::name; }
};

} // namespace clem
#endif
