/**
 * The eddies of the stirring: their range of sizes on a level, and the
 * random draws of one stirring event.
 *
 * The size l of the eddy of an event follows the inertial-range PDF
 *
 *   f(l) = (5/3) l^(-8/3) / ( l_min^(-5/3) - Delta^(-5/3) ),   l_min <= l <= Delta
 *
 * drawn by inverting its CDF, F(l) = (l_min^(-5/3) - l^(-5/3)) / (l_min^(-5/3) - Delta^(-5/3)):
 *
 *   l = ( l_min^(-5/3) - u (l_min^(-5/3) - Delta^(-5/3)) )^(-3/5),   u uniform on [0, 1)
 *
 * Every element is taken as l_particle = L / n_lem long (L the cell size
 * along x), rather than by its own extent in the mass coordinate. The eddy
 * then covers l / l_particle elements, rounded to the nearest multiple of 3
 * (a triplet map lays down three compressed copies of the segment, so the
 * segment splits into three equal thirds): 5.5 -> 6, 7 -> 6, 8 -> 9. The
 * smallest map is setup::stirring_min_elements = 6 (a map of 3 is the
 * identity), and f(l) is truncated there, l_min >= 6 l_particle, so every draw
 * gives a map that acts. Its place on the line is uniform:
 * the first element of the segment is any of the n_elem with the same
 * probability, the line being periodic.
 *
 * Functions:
 *
 * LevelEddyRange  : the sizes of a level: L, l_particle, Delta and the resolution limit
 *
 * SampleEddySize  : l from f(l), given a uniform u
 *
 * EddyElements    : the number of elements an eddy of size l covers
 *
 * SampleEddyStart : the first element of the segment, given a uniform u
 */
#ifndef CLEM_EDDY_SAMPLING_H
#define CLEM_EDDY_SAMPLING_H

#include <cmath>

#include <AMReX_Geometry.H>

#include "setup.H"

namespace clem {

// the sizes of the eddies of the lines of one level
struct EddyRange
{
  amrex::Real line_length = 0.0; // L, the cell size along x
  amrex::Real l_particle  = 0.0; // the length of an element, L / n_lem
  amrex::Real delta       = 0.0; // the largest eddy, min(filter size, L)
  amrex::Real l_resolved  = 0.0; // the smallest eddy the line holds,
                                 // stirring_min_elements l_particle

  // the smallest eddy of a cell of Kolmogorov scale eta
  amrex::Real Smallest(const amrex::Real eta) const { return amrex::max(eta, l_resolved); }

  // does a cell of Kolmogorov scale eta have an inertial range on its line
  bool Stirs(const amrex::Real eta) const { return eta > 0.0 && Smallest(eta) < delta; }
};

// LevelEddyRange: the eddy sizes of the level of geom. The filter size is
// (dx dy dz)^(1/3), as in the WALE closure
inline EddyRange
LevelEddyRange(const amrex::Geometry& geom)
{
  const auto cell_size = geom.CellSizeArray();
  EddyRange range;
  range.line_length = cell_size[0];
  range.l_particle  = range.line_length / static_cast<amrex::Real>(setup::n_lem);
  const amrex::Real filter_size =
    std::pow(AMREX_D_TERM(cell_size[0], *cell_size[1], *cell_size[2]), 1.0 / AMREX_SPACEDIM);
  range.delta      = amrex::min(filter_size, range.line_length);
  range.l_resolved = static_cast<amrex::Real>(setup::stirring_min_elements) * range.l_particle;
  return range;
}

// SampleEddySize: the size of an eddy drawn from f(l) on [l_min, delta], by
// the inverse of its CDF at u in [0, 1)
inline amrex::Real
SampleEddySize(const amrex::Real l_min, const amrex::Real delta, const amrex::Real u)
{
  const amrex::Real a = std::pow(l_min, -5.0 / 3.0);
  const amrex::Real b = std::pow(delta, -5.0 / 3.0);
  return std::pow(a - u * (a - b), -3.0 / 5.0);
}

// EddyElements: the elements an eddy of size l covers, l / l_particle rounded
// to the nearest multiple of 3, kept within [stirring_min_elements, the
// largest multiple of 3 that fits in the n_elem elements of the line]. 0 if
// the line is too short for any triplet map
inline int
EddyElements(const amrex::Real l, const amrex::Real l_particle, const int n_elem)
{
  const int n_max = 3 * (n_elem / 3);
  if (n_max < setup::stirring_min_elements) {
    return 0;
  }
  const amrex::Real n_exact = l / l_particle;
  const int n_seg           = 3 * static_cast<int>(std::lround(n_exact / 3.0));
  return amrex::min(amrex::max(n_seg, setup::stirring_min_elements), n_max);
}

// SampleEddyStart: the first element of the segment of the eddy, uniform on
// the n_elem elements of the periodic line, at u in [0, 1)
inline int
SampleEddyStart(const int n_elem, const amrex::Real u)
{
  const int start = static_cast<int>(u * static_cast<amrex::Real>(n_elem));
  return amrex::min(amrex::max(start, 0), n_elem - 1);
}

} // namespace clem
#endif
