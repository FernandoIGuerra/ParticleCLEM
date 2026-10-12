/**
 * This is a derived class from Reindexer.
 *
 * AlternatingReindexer: omega flips every step of the level,
 *
 *   ... o R o Sigma^+ o R o Sigma^- ,
 *
 * Sigma^- on the even steps and Sigma^+ on the odd ones, each over the full
 * transfer of its step and each followed by its regrid. The donor offset and
 * the position of a block among the arrivals then average to the centre of
 * their range over two steps.
 */
#ifndef CLEM_ALTERNATING_REINDEXER_H
#define CLEM_ALTERNATING_REINDEXER_H

#include "reindexer.H"

namespace clem {

class AlternatingReindexer : public Reindexer<AlternatingReindexer>
{
  friend Reindexer;

protected:
  // O^omega on N^+_c
  void SortFluxImplementation(int* face, const amrex::Real* mu, const int m, const int step) const
  {
    detail::OrderTransfers(face, mu, m, step % 2 != 0);
  }

  // O^omega on N^-_c
  void SortIncomingParticleImplementation(int* face, const amrex::Real* mu, const int m,
                                          const int step) const
  {
    detail::OrderTransfers(face, mu, m, step % 2 != 0);
  }

public:
  static constexpr const char* name = "Sigma^- / Sigma^+ alternating every step";

  AlternatingReindexer()                     = default;
  ~AlternatingReindexer()                    = default;
  AlternatingReindexer(AlternatingReindexer const&)            = delete;
  AlternatingReindexer& operator=(AlternatingReindexer const&) = delete;
  AlternatingReindexer(AlternatingReindexer&&)                 = delete;
  AlternatingReindexer& operator=(AlternatingReindexer&&)      = delete;
};

} // namespace clem
#endif
