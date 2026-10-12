/**
 * Here we create the base class of the ordering policy of the splicing.
 *
 * Reindexer<Derived>: the base of every policy (CRTP, no virtual calls). A
 * cell has outflux faces N^+_c and influx faces N^-_c, each with its transfer
 * mass mu_f of the step. The policy is the order O^omega in which they are
 * taken:
 *
 *   SortFlux              O^omega on N^+_c: the schedule in which a donor
 *                         serves its outflux faces, each from the head of
 *                         what is left of its line
 *   SortIncomingParticle  O^omega on N^-_c: the rank in which a receptor
 *                         appends the arriving blocks at its tail
 *
 * The SAME omega on both sides is what makes the order in which material
 * leaves one line the order in which it enters the next.
 *
 * step is the step count of the level, for the policies that change with it.
 *
 * The policy is chosen at compile time in setup.H.
 */
#ifndef CLEM_REINDEXER_H
#define CLEM_REINDEXER_H

#include <algorithm>
#include <cmath>

#include <AMReX_REAL.H>

namespace clem {

namespace detail {
// sort face[0..m) by the SIZE of the transfer, |mu|, largest first if
// descending; ties by face index, so the order is a strict total order.
// The size, not the signed value: seen from a receptor the transfers are
// negative, and Sigma^+ must still append the largest block first. With
// arrivals of 0.5 from A and 1 from B, the order is B then A whether they
// are given as (0.5, 1) or as (-0.5, -1)
inline void
OrderTransfers(int* face, const amrex::Real* mu, const int m, const bool descending)
{
  std::sort(face, face + m, [=](const int a, const int b) {
    const amrex::Real size_a = std::abs(mu[a]);
    const amrex::Real size_b = std::abs(mu[b]);
    if (size_a != size_b) {
      return descending ? (size_a > size_b) : (size_a < size_b);
    }
    return a < b;
  });
}
} // namespace detail

template <typename DerivedReindexer>
class Reindexer
{
  friend DerivedReindexer;

  // only a derived policy can build the base
  Reindexer() = default;

  const DerivedReindexer& derived() const { return static_cast<const DerivedReindexer&>(*this); }

public:
  ~Reindexer()                           = default;
  Reindexer(Reindexer const&)            = delete;
  Reindexer& operator=(Reindexer const&) = delete;
  Reindexer(Reindexer&&)                 = delete;
  Reindexer& operator=(Reindexer&&)      = delete;

  // reorder face[0..m) into the schedule of the donor; mu[f] is the transfer
  // mass of face f
  void SortFlux(int* face, const amrex::Real* mu, const int m, const int step) const
  {
    derived().SortFluxImplementation(face, mu, m, step);
  }

  // reorder face[0..m) into the rank in which the receptor appends the
  // blocks; mu[f] is the transfer of face f, of either sign: only its size
  // is used
  void SortIncomingParticle(int* face, const amrex::Real* mu, const int m, const int step) const
  {
    derived().SortIncomingParticleImplementation(face, mu, m, step);
  }

  static constexpr const char* Name() { return DerivedReindexer::name; }
};

} // namespace clem
#endif
