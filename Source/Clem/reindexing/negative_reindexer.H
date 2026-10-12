/**
 * This is a derived class from Reindexer.
 *
 * NegativeReindexer: Sigma^-, the smallest transfer first, on every step.
 */
#ifndef CLEM_NEGATIVE_REINDEXER_H
#define CLEM_NEGATIVE_REINDEXER_H

#include "reindexer.H"

namespace clem {

class NegativeReindexer : public Reindexer<NegativeReindexer>
{
  friend Reindexer;

protected:
  // O^omega on N^+_c
  void SortFluxImplementation(int* face, const amrex::Real* mu, const int m, const int /*step*/) const
  {
    detail::OrderTransfers(face, mu, m, false);
  }

  // O^omega on N^-_c
  void SortIncomingParticleImplementation(int* face, const amrex::Real* mu, const int m,
                                          const int /*step*/) const
  {
    detail::OrderTransfers(face, mu, m, false);
  }

public:
  static constexpr const char* name = "Sigma^- (smallest transfer first)";

  NegativeReindexer()                     = default;
  ~NegativeReindexer()                    = default;
  NegativeReindexer(NegativeReindexer const&)            = delete;
  NegativeReindexer& operator=(NegativeReindexer const&) = delete;
  NegativeReindexer(NegativeReindexer&&)                 = delete;
  NegativeReindexer& operator=(NegativeReindexer&&)      = delete;
};

} // namespace clem
#endif
