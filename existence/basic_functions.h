/* Blah blah blah. Javi */

#ifndef BASIC_FUNCTIONS_H_INCLUDED
#define BASIC_FUNCTIONS_H_INCLUDED

#include <string.h>
#include "arb_calc.h"
#include "arb_poly.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "arb_hypgeom.h"
#include "flint/profiler.h"

#define TANH_MAX 32

typedef struct {
  int N0,N1,N2;
  arb_t p0,p1;
  arb_ptr bi;
  arb_ptr a0k;
  arb_ptr a1k;
} u_params_struct;

typedef struct {
  u_params_struct up;
  arb_t x;
  arb_ptr u;
  int order;
} Integration_params_struct;


typedef int (*whitham_calc_func)(arb_t out, const arb_t inp, void * param, slong prec);
typedef int (*explicit_integral)(arb_t out, void * param, slong prec);
typedef int (*explicit_singularity)(arb_t out, const arb_t delta, void * param, slong prec);

// Basic functions related to the special functions:

// Clausen functions
void arb_clausen_sin_direct(arb_t out, const arb_t s, const arb_t x, slong prec);
void arb_clausen_cos(arb_t out, const arb_t s, const arb_t x, slong prec);
void arb_clausen_sin(arb_t out, const arb_t s, const arb_t x, slong prec);

// Cosintegral function
// Due to the suboptimal performance of the arb function for big real intervals,
// we combine it with a Taylor expansion at infinity
// 1/x*(sin(x)-1/x*(cos(x)+1)) <= Ci(x) <= 1/x*(sin(x)-1/x*(cos(x)-1))
void arb_cosintegral(arb_t out, const arb_t x, slong prec);


// TODO(Javi): Maybe this is not needed!!
// Fresnel-sine/cosine function
// Due to the suboptimal performance of the arb function for big real intervals,
// we combine it with a Taylor expansion at infinity
// ... ?
void arb_fresnel_sin_cos(arb_t out1, arb_t out2, const arb_t x, slong prec);


// Basic functions related to tail estimates. Main estimate.
// Computes error bounds on \sum_{n>N0} (\sqrt{tanh(nx)}-1)*n^{exponent}. See Appendix B.
void compute_tail_error(arb_ptr out, const arb_t exponent, const int N0, slong prec);



// Computation of the functions related to u0

// Populates a u_params_struct with the right parameters for u0
void populate_parameters(u_params_struct* up, slong prec);

// Populates the spatial part of an Integration_params_struct with x and the values of u(x),u'(x),... up to order 'order'
void populate_Integration_parameters(Integration_params_struct* ips, const arb_t x, int order, slong prec);



int compute_u0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_du0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_ddu0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dddu0(arb_ptr out, const arb_t inp, void* params, slong prec);

int compute_Lu0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dLu0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_ddLu0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dddLu0(arb_ptr out, const arb_t inp, void* params, slong prec);

// Computation of the functions related to the kernels: its derivatives and integrals

void IIKp(arb_t out, const arb_t x, slong prec);
void IKp(arb_t out, const arb_t x, slong prec);
void Kp(arb_t out, const arb_t x, slong prec);
void dKp(arb_t out, const arb_t x, slong prec);
void ddKp(arb_t out, const arb_t x, slong prec);



// For debugging purposes
int zero(arb_ptr out, const arb_t delta, void* params, slong prec);




#endif
