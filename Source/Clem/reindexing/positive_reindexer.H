/**
 * This is a derived class from Reindexer.
 *
 * PositiveReindexer: Sigma^+, the largest transfer first, on every step.
 */
#ifndef CLEM_POSITIVE_REINDEXER_H
#define CLEM_POSITIVE_REINDEXER_H

#include "reindexer.H"

namespace clem {

class PositiveReindexer : public Reindexer<PositiveReindexer>
{
  friend Reindexer;

protected:
  // O^omega on N^+_c
  void SortFluxImplementation(int* face, const amrex::Real* mu, const int m, const int /*step*/) const
  {
    detail::OrderTransfers(face, mu, m, true);
  }

  // O^omega on N^-_c
  void SortIncomingParticleImplementation(int* face, const amrex::Real* mu, const int m,
                                          const int /*step*/) const
  {
    detail::OrderTransfers(face, mu, m, true);
  }

public:
  static constexpr const char* name = "Sigma^+ (largest transfer first)";

  PositiveReindexer()                     = default;
  ~PositiveReindexer()                    = default;
  PositiveReindexer(PositiveReindexer const&)            = delete;
  PositiveReindexer& operator=(PositiveReindexer const&) = delete;
  PositiveReindexer(PositiveReindexer&&)                 = delete;
  PositiveReindexer& operator=(PositiveReindexer&&)      = delete;
};

} // namespace clem
#endif
