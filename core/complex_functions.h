/* complex_functions.h: pointwise evaluators built on u_0 and the kernel.
 *
 * Every routine stores a ball in out and returns 0; the return value is not
 * used by any program (see basic_functions.h). The argument inp is the point
 * x unless stated otherwise.
 */
#ifndef COMPLEX_FUNCTIONS_H_INCLUDED
#define COMPLEX_FUNCTIONS_H_INCLUDED

#include "basic_functions.h"

/* The coefficient functions of the lemma on the coefficient bounds, whose
   suprema over (0, pi] are alpha_1, alpha_2 and alpha_f (derivative_bounds.c
   bounds them on the compact cells; derivative_cusp.c handles the cusp). */
int compute_alpha1(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_alpha2(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_alphaf(arb_ptr out, const arb_t inp, void* params, slong prec);

/* The curvature function x u_0''(x). */
int compute_convexity(arb_ptr out, const arb_t inp, void* params, slong prec);

/* The defect F_0 of the approximate solution and its derivative. */
int compute_F0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dF0(arb_ptr out, const arb_t inp, void* params, slong prec);

/* The explicit integral of the T_1 integrand over the far region y > x
   (Appendix B); params points to an Integration_params_struct holding x and
   the values of u_0 and its derivatives at x. */
int compute_big_y_T1(arb_ptr out, void* params, slong prec);

#endif
