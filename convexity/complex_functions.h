/* Blah blah blah. Javi */

#ifndef COMPLEX_FUNCTIONS_H_INCLUDED
#define COMPLEX_FUNCTIONS_H_INCLUDED

#include "basic_functions.h"

// TODO(Javi): Mirar bien todos los return 0 vs return 1. Estan descompensados y en unos sitios 0 es bueno, en otros 1 es bueno.!!!!
// Otros son void, etc. Menudo desastre.


int compute_alpha0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_alpha1(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_alpha2(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_alpha2bar(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_alphaf(arb_ptr out, const arb_t inp, void* params, slong prec);

int compute_convexity(arb_ptr out, const arb_t inp, void* params, slong prec);

// inp = w, x comes from params
void compute_K1(arb_ptr out, const arb_t inp, void* params, slong prec);
void compute_K20(arb_ptr out, const arb_t inp, void* params, slong prec);
void compute_K2bar0(arb_ptr out, const arb_t inp, void* params, slong prec);

void compute_dyK20(arb_ptr out, const arb_t inp, void* params, slong prec);
void compute_dyK2bar0(arb_ptr out, const arb_t inp, void* params, slong prec);
void compute_ddyK20(arb_ptr out, const arb_t inp, void* params, slong prec);
void compute_ddyK2bar0(arb_ptr out, const arb_t inp, void* params, slong prec);


int compute_F0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dF0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_ddF0(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dddF0(arb_ptr out, const arb_t inp, void* params, slong prec);

int compute_F1(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dF1(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_ddF1(arb_ptr out, const arb_t inp, void* params, slong prec);


// TODO(Javi): These two may be optimized by implementing them via F0 and derivatives instead of F1
int compute_F2(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_dF2(arb_ptr out, const arb_t inp, void* params, slong prec);

// Far regions
// We calculate Ti in the region y > x
int compute_big_y_T0(arb_ptr out, void* params, slong prec);
int compute_big_y_T1(arb_ptr out, void* params, slong prec);
int compute_big_y_T2(arb_ptr out, void* params, slong prec);


// Integrands

int integrand_T0(arb_ptr out, const arb_t inp, void* params, slong prec);
//int integrand_dT0(arb_ptr out, const arb_t inp, void* params, slong prec); //TODO(Javi): Do we need this one?
int integrand_ddT0(arb_ptr out, const arb_t inp, void* params, slong prec);

int integrand_T1(arb_ptr out, const arb_t inp, void* params, slong prec);
int integrand_ddT1(arb_ptr out, const arb_t inp, void* params, slong prec);

int integrand_T2(arb_ptr out, const arb_t inp, void* params, slong prec);
int integrand_ddT2(arb_ptr out, const arb_t inp, void* params, slong prec);

// Singularities

int singularity_T0_0(arb_ptr out, const arb_t delta, void* params, slong prec);
int singularity_T0_1(arb_ptr out, const arb_t delta, void* params, slong prec);
int singularity_T1_1(arb_ptr out, const arb_t delta, void* params, slong prec);
int singularity_T2_0(arb_ptr out, const arb_t delta, void* params, slong prec);
int singularity_T2_1(arb_ptr out, const arb_t delta, void* params, slong prec);

// Auxiliary functions for the singularities
int compute_C_Delta(arb_ptr out, const arb_t delta, slong prec);
int compute_c_delta(arb_ptr out, const arb_t delta, slong prec);




#endif
