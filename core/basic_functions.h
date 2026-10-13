/* basic_functions.h: shared evaluators of the computer-assisted proof.
 *
 * Declares the special functions (Clausen functions, through the cached
 * series of clausen_cache.c), the Fourier-tail bound of Appendix B, the
 * approximate solution u_0 of equation (u0) with its derivatives, the
 * operator L applied to u_0, and the primitives and derivatives of the
 * kernel. All routines use FLINT ball arithmetic at the precision passed as
 * prec; every program passes CAP_PREC (precision.h).
 *
 * The compute_* routines share the signature of a callback: they store an
 * enclosure in out and return an int. No program in this package uses the
 * return value; the index routines of constant_functions.c return 0 for an
 * index outside the expansion and 1 otherwise, the other routines return 0.
 */
#ifndef BASIC_FUNCTIONS_H_INCLUDED
#define BASIC_FUNCTIONS_H_INCLUDED

#include <string.h>
#include "arb_calc.h"
#include "arb_poly.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "arb_hypgeom.h"

/* Fourier sums with the weights sqrt(tanh n) are truncated after n = TANH_MAX;
   the tails are bounded from n = TANH_MAX + 1 by compute_tail_error. */
#define TANH_MAX 32

/* The approximate solution u_0 of equation (u0): the exponents p0, p1 and the
   coefficients a0k (N0 terms), a1k (N1 terms) and bi (N2 terms), exact
   decimals set by populate_parameters (inputs/approximation.json). */
typedef struct {
  int N0,N1,N2;
  arb_t p0,p1;
  arb_ptr bi;
  arb_ptr a0k;
  arb_ptr a1k;
} u_params_struct;

/* Parameters of the kernel integrands: the approximation, a point x and the
   values of u_0 and its derivatives at x up to the given order. */
typedef struct {
  u_params_struct up;
  arb_t x;
  arb_ptr u;
  int order;
} Integration_params_struct;

/* Clausen functions sum_{n>=1} sin(nx)/n^s and sum_{n>=1} cos(nx)/n^s for
   real orders s. clausen_cached (clausen_cache.c) evaluates the series of
   Appendix B for non-integer s and returns 1 when it has done so; the
   direct evaluators handle the remaining cases. */
int clausen_cached(arb_t out, const arb_t s, const arb_t x, int sine, slong prec);
void arb_clausen_sin_direct(arb_t out, const arb_t s, const arb_t x, slong prec);
void arb_clausen_cos(arb_t out, const arb_t s, const arb_t x, slong prec);
void arb_clausen_sin(arb_t out, const arb_t s, const arb_t x, slong prec);

/* Fourier-tail bound of Appendix B: an upper bound for
   sum_{n>=N0} n^exponent (1 - sqrt(tanh n)), valid for every value of the
   exponent in its ball. The lemma requires N0 >= 10 and the positivity of
   2 - max(exponent, 0) log(N0)/N0, which is asserted. */
void compute_tail_error(arb_ptr out, const arb_t exponent, const int N0, slong prec);

/* Sets the exact coefficients of u_0. */
void populate_parameters(u_params_struct* up, slong prec);

/* u_0 and its first three derivatives at inp; L u_0 and its first two
   derivatives at inp. */
int compute_u0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_du0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_ddu0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_Lu0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dLu0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_ddLu0(arb_ptr out, const arb_t inp, void* params, slong prec);

/* The kernel K of the paper through its Fourier series: Kp(x) = K(x), its
   first and second primitives IKp and IIKp vanishing at 0 (IIKp with
   vanishing first derivative as well), and its derivatives dKp and ddKp.
   Each is a Clausen function of order 1/2, 3/2, 5/2, -1/2 or -3/2 plus the
   finite correction with the weights sqrt(tanh n) - 1 and the tail bound. */
void IIKp(arb_t out, const arb_t x, slong prec);
void IKp(arb_t out, const arb_t x, slong prec);
void Kp(arb_t out, const arb_t x, slong prec);
void dKp(arb_t out, const arb_t x, slong prec);
void ddKp(arb_t out, const arb_t x, slong prec);

#endif
