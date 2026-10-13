/* Blah blah blah. Javi */

#ifndef METHODS_H_INCLUDED
#define METHODS_H_INCLUDED

#include <string.h>
#include "arb_calc.h"
#include "arb_poly.h"
#include "acb_hypgeom.h"
#include "basic_functions.h"
#include "complex_functions.h"

// TODO(Javi): Maybe put a flag that says VERBOSE and print or not depending on this flag.


// Checks if the function "func" is smaller than the bound "bound".
// We subdivide until the minimum_width "min_width".
// We use C^0 information
int compute_bound_C0(const arb_t min_width, const arb_t bound, whitham_calc_func func, void * param, const arb_t inp, slong prec);



// Checks if the function "func" is smaller than the bound "bound".
// We subdivide until the minimum_width "min_width".
// We use C^1 information
int compute_bound_C1(const arb_t min_width, const arb_t bound, whitham_calc_func func, whitham_calc_func deriv, void * param, const arb_t inp, slong prec);




////////////////////////////// OLD STUFF

// Input: an IntegrationResult IR with all the system parameters (passed by reference),
//        a function integrand_1D_fast which is interval-based,
//        and a function integrand_1D_slow which is itaylor-based:
//        integrand_1D_fast and integrand_1D_slow represent the integrand.
// Output: the result and the error will be stored in IR
// The function computes a Gauss-Legendre quadrature of order 2 between IR.p.Left and IR.p.Right.
// To optimize the code, we only call integrand_1D_slow whenever we want to
// get information about the derivatives of the integrand, otherwise we
// call integrand_1D_fast.
// IR.flag is set to 1 if the integrands raise a flag.
// Note that this function returns the value of the integral divided by the length

//typedef int (*whitham_calc_func)(arb_t out, const arb_t inp, void * param, slong prec);
// Integrand_Function is the main term
// Integrand_Derivative is the error term
// returns 0 if OK, 1 if error
// Do Midpoint if possible, otherwise C^0
int MidPointQuadratureAbs(arb_t out, const arb_t inp_left, const arb_t inp_right, void* param, 
                       whitham_calc_func Integrand_Function, whitham_calc_func Integrand_Derivative, slong prec);



// 1D callers

// Input: a vector of ParameterSets called vp, containing the set of regions over which we
//        will integrate,
//        an integration method (in this case only new_GaussLegendre_2_1D is available),
//        a function integrand_1D_fast which is interval-based,
//        a function integrand_1D_slow which is itaylor-based:
//        integrand_1D_fast and integrand_1D_slow represent the integrand,
//        and an integer flag signaling possible errors.
// Output: an interval containing the sum of the integrals over all members of vp.
// This function performs an adaptive integration, integrating first over
// the intervals that yield biggest enclosures in absolute terms. The intervals
// are stored in a priority_queue in a sorted way.
// Integration regions are accepted or rejected based on absolute and relative
// tolerances (with respect to the length of the domain), encoded in each ParameterSet.
// If a region is rejected, it is split by the midpoint of the interval and the two
// resulting intervals are pushed again in the priority_queue.
// We repeat this procedure up to MAX_ELEMENTS_EVALUATED times. If we reach the limit,
// we sum the result anyway and set the second rightmost bit of 'flag' to 1.

// TODO(Javi): De momento dejalo asi, mas adelante quiza lo puedes poner como una Priority Queue.
// Quito el "metodo de integracion" porque no lo cambiaremos.
int integrate(arb_t out, const arb_t inp_left, const arb_t inp_right, void* param, 
              whitham_calc_func Integrand_Function, whitham_calc_func Integrand_Derivative, slong prec,
              const arb_t max_width_abs, const arb_t max_width_rel, int error_mask);



// Checks if the Linfinity norm is smaller than the bound "bound".
// We subdivide (in x) until the minimum_width "min_width".
int bound_operator(const arb_t inp_left, const arb_t inp_right, void* param, 
                   whitham_calc_func Integrand_Function, whitham_calc_func Integrand_Derivative,
                   explicit_integral BigY_Function, explicit_singularity SingularityZero_Function,
                   explicit_singularity SingularityOne_Function, const arb_t delta0, const arb_t delta1,
                   slong prec, const arb_t min_width_x, const arb_t max_width_abs, const arb_t max_width_rel, 
                   int error_mask, const arb_t bound);


#endif
