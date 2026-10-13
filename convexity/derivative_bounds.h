#ifndef DERIVATIVE_BOUNDS_H
#define DERIVATIVE_BOUNDS_H
#include "complex_functions.h"
/* Correct derivatives of R=V/(2U)-U/2: F1=R', F2=x R''. */
void derivative_defect_jets(arb_ptr out, const arb_t x, u_params_struct *up, int order, slong p);
void derivative_cusp_errors(arb_ptr out, const arb_t eps, u_params_struct *up, slong prec);
#endif
