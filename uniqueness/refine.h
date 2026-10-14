/* refine.h: exact integer refinement of the a priori enclosures (refine.c).
 *
 * V is a 64-bit lattice value, I a signed 128-bit accumulator. A Sum holds
 * prefix sums of a nonnegative row against the node values; a Row gathers
 * the three rows of one node with its clipped minimizers. advance() performs
 * one coupled refinement step on all nodes; the refine_* routines are its
 * building blocks, exported so that test_refine.c can check them exhaustively.
 */
#ifndef REFINE_H_INCLUDED
#define REFINE_H_INCLUDED
#include <stdint.h>
typedef __int128 I;
typedef int64_t V;
typedef struct {
  int n;
  const V *v;
  V *s;
  I *sv;
} Sum;
typedef struct {
  Sum a, b, m, c;
  V pu, tu;
  I err, base, outer;
  int difference;
} Row;
I refine_square(V q);
Sum refine_make(const V *a, const V *v, int n);
void refine_clear(Sum *s);
V refine_minp(const Row *r, V lo, V hi);
V refine_mint(const Row *r, V lo, V hi);
I refine_inner(const Row *r, V t);
int advance(int n, const V *fl, const V *fh, const V *gl, const V *gh, const V *dpl, const V *dnh,
            const int *zi, V *L, V *U, V l0, V h0, int differences);
#endif
