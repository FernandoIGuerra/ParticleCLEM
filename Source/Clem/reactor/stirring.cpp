/**
 * stirring.cpp: the stirring of one line.
 *
 *   DoStirring  a stirring event: the triplet map of one eddy on the line
 *
 * The map of an eddy of l elements (a multiple of 3) starting at element j
 * (Kerstein 1991) rearranges the segment into three compressed copies of
 * itself, the middle one reversed:
 *
 *   j, j+3, j+6, ..., j+l-3,   j+l-2, j+l-5, ..., j+4, j+1,   j+2, j+5, ..., j+l-1
 *
 * that is, place i of the segment takes the element at offset
 *
 *   3i             first third    i = 0 .. l/3 - 1
 *   l - 2 - 3i     middle third   (reversed)
 *   2 + 3i         last third
 *
 * from j. The line is periodic: the offsets wrap around its n_elem elements.
 * It is a permutation: every element keeps its mass and its state, so the
 * map conserves everything exactly. The eddy (l, j) is drawn outside
 * (stirring/eddy_sampling.H).
 */
#include "line_operators.H"

namespace clem {

namespace {

// the segment of `field` (stride values per element) rearranged: place o
// (element (start + o) mod n) takes element (start + source[o]) mod n,
// through `buffer`
template <typename T>
void
PermuteSegment(amrex::Vector<T>& field, const int stride, const amrex::Vector<int>& source,
               const int start, const int n_eddy, const int n_elem, amrex::Vector<T>& buffer)
{
  for (int o = 0; o < n_eddy; ++o) {
    const int from = (start + source[o]) % n_elem;
    for (int c = 0; c < stride; ++c) {
      buffer[o * stride + c] = field[from * stride + c];
    }
  }
  for (int o = 0; o < n_eddy; ++o) {
    const int to = (start + o) % n_elem;
    for (int c = 0; c < stride; ++c) {
      field[to * stride + c] = buffer[o * stride + c];
    }
  }
}

} // namespace

void
DoStirring(LineScratch& line, const int start, const int n_eddy)
{
  const int n_elem = line.n_elem;
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_eddy % 3 == 0 && n_eddy >= setup::stirring_min_elements &&
                                     n_eddy <= n_elem && start >= 0 && start < n_elem,
                                   "DoStirring: the eddy does not fit the line");

  // the source offset of every place of the segment
  const int third = n_eddy / 3;
  auto& source    = line.stir_source;
  for (int i = 0; i < third; ++i) {
    source[i]             = 3 * i;              // j, j+3, ..., j+l-3
    source[third + i]     = n_eddy - 2 - 3 * i; // j+l-2, j+l-5, ..., j+1
    source[2 * third + i] = 2 + 3 * i;          // j+2, j+5, ..., j+l-1
  }

  // the stored state of the elements, and which particle each one is; the
  // derived state (X, Cp, transport) follows at the next UpdateProperties
  PermuteSegment(line.particle_index, 1, source, start, n_eddy, n_elem, line.stir_int);
  PermuteSegment(line.mass, 1, source, start, n_eddy, n_elem, line.stir_real);
  PermuteSegment(line.rho, 1, source, start, n_eddy, n_elem, line.stir_real);
  PermuteSegment(line.T, 1, source, start, n_eddy, n_elem, line.stir_real);
  PermuteSegment(line.eint, 1, source, start, n_eddy, n_elem, line.stir_real);
  PermuteSegment(line.press, 1, source, start, n_eddy, n_elem, line.stir_real);
  PermuteSegment(line.Y, NUM_SPECIES, source, start, n_eddy, n_elem, line.stir_real);
}

} // namespace clem
