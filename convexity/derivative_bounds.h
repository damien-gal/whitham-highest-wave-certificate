/* derivative_bounds.h: derivative jets shared by the convexity programs.
 * derivative_defect_jets() encloses the derivatives of R = V/(2U) - U/2, so
 * that F_1 = R' and F_2 = x R''; derivative_cusp_errors() gives the cusp
 * remainder bounds of derivative_cusp.c.
 */
#ifndef DERIVATIVE_BOUNDS_H
#define DERIVATIVE_BOUNDS_H
#include "complex_functions.h"
/* Correct derivatives of R=V/(2U)-U/2: F1=R', F2=x R''. */
void derivative_defect_jets(arb_ptr out, const arb_t x, u_params_struct *up, int order, slong p);
void derivative_cusp_errors(arb_ptr out, const arb_t eps, u_params_struct *up, slong prec);
#endif
