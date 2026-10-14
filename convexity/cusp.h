/* cusp.h: the cusp expansions shared by the convexity programs (cusp.c).
 *
 * On 0 < x <= eps the scaled profile U = u_0/sqrt(x) and the residual of the
 * equation are finite sums of powers x^(h/2 + k p_0 + j p_1) with explicit
 * remainders. cusp_init() fixes the approximation and builds the coefficient
 * tables once; the other routines evaluate exponents, falling factorials,
 * the complete tanh sums, the Clausen tails and the cusp remainders, and
 * derivative_cusp_errors() and derivative_profile_jets() are the interfaces
 * used by the local estimates and the operator programs.
 */
#ifndef CUSP_H_INCLUDED
#define CUSP_H_INCLUDED
#include "constant_functions.h"

#define CUSP_M 7   /* even terms kept in the regular parts */
#define CUSP_N 128 /* terms of the finite tanh sums */

/* One term c x^(h/2 + k p_0 + j p_1) of an expansion. */
typedef struct {
  int h, k, j;
  arb_t c;
} cusp_term;

extern u_params_struct cusp_up;          /* the approximation fixed by cusp_init */
extern cusp_term cusp_us[24], cusp_ss[512]; /* terms of U and of the residual */
extern int cusp_nu, cusp_ns;             /* their numbers */

void cusp_init(const u_params_struct *params);
void cusp_exponent(arb_t a, int h, int k, int j);
void cusp_falling(arb_t a, const arb_t e, int d);
void cusp_smooth(arb_t out, const arb_t e);
void cusp_cltail(arb_t out, const arb_t s, const arb_t x, int d);
void cusp_remainder(arb_t out, const arb_t x, int d, int v);
void derivative_cusp_errors(arb_ptr out, const arb_t eps, u_params_struct *params, slong prec);
int derivative_profile_jets(arb_ptr out, const arb_t X, u_params_struct *params, slong prec);

#endif
