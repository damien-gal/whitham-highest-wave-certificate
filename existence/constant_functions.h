/* Blah blah blah. Javi */

#ifndef CONSTANT_FUNCTIONS_H_INCLUDED
#define CONSTANT_FUNCTIONS_H_INCLUDED

#include "basic_functions.h"


int compute_a_jk_0(arb_ptr out, int j, int k, void* params, slong prec);
int compute_a_jk_1(arb_ptr out, int j, int k, void* params, slong prec);
int compute_a_jk_2(arb_ptr out, int j, int k, void* params, slong prec);

int compute_a_m_0(arb_ptr out, int m, void* params, slong prec);
int compute_a_m_1(arb_ptr out, int m, void* params, slong prec);
int compute_a_m_2(arb_ptr out, int m, void* params, slong prec);

int compute_A_jk_0(arb_ptr out, int j, int k, void* params, slong prec);

int compute_A_m_0(arb_ptr out, int m, void* params, slong prec);
int compute_A_m_1(arb_ptr out, int m, void* params, slong prec);
int compute_A_m_2(arb_ptr out, int m, void* params, slong prec);

int compute_b_0_0(arb_ptr out, void* params, slong prec);
int compute_b_0_1(arb_ptr out, void* params, slong prec);
int compute_b_0_2(arb_ptr out, void* params, slong prec);
int compute_b_1_0(arb_ptr out, void* params, slong prec);
int compute_b_1_1(arb_ptr out, void* params, slong prec);
int compute_b_1_2(arb_ptr out, void* params, slong prec);

int compute_E_x_Cz(arb_ptr out, const arb_t inp, const arb_t z, void* params, slong prec);
int compute_E_x_Sz(arb_ptr out, const arb_t inp, const arb_t z, void* params, slong prec);

int compute_E_x_u0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_E_x_du0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_E_x_ddu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);

int compute_E_x_Lu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_E_x_dLu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_E_x_ddLu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);

int compute_E_x_f(arb_ptr out, const arb_t inp, int M, void* params, slong prec);

int compute_E_x_F0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_E_x_F1(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_E_x_F2(arb_ptr out, const arb_t inp, int M, void* params, slong prec);

int compute_c_hat_plus(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_c_hat_minus(arb_ptr out, const arb_t inp, int M, void* params, slong prec);

int compute_c_u0_hat_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_u0_hat_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

int compute_c_p0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

int compute_cp_p0_zero(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_cp_p0_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_cp_p0_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

int compute_cpp_p0_zero(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_cpp_p0_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_cpp_p0_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);


int compute_c_F0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_F1(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_F2(arb_ptr out, const arb_t eps, int M, void* params, slong prec);


int compute_bound_alpha_eps_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_bound_alpha1(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_bound_alpha2(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_bound_alpha2bar(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_bound_alphaf(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

int compute_bound_delta_eps_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_bound_delta_eps_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_bound_delta_eps_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec);


// TODO(Javi): inp vs eps???
int compute_cT0_1(arb_ptr out, const arb_t inp, void* params, slong prec);
int compute_cT0_2(arb_ptr out, const arb_t inp, void* params, slong prec);

int compute_cT1_1(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_cT1_2(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_cT1_3(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_cT1_4(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_cp_eps(arb_ptr out, const arb_t inp, int M, void* params, slong prec);


int compute_c_p0_SE_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_p0_SE_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

int compute_c_p0_T2_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_p0_T2_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_p0_T2_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

int compute_cT2_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_cT2_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_cT2_2_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_cT2_3(arb_ptr out, const arb_t eps, int M, void* params, slong prec);




// (XXX)
int compute_validation_du0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

// (5.10) in the old version
int compute_validation_T0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

// (6.24) in the paper. 
// We prove (c_{ep,p0} + c'_{ep,p0,-}(c_{T_1}^{1} + c_{T_1}^{2}))\ep^{p0-1/2} + c_{T_1}^{3} \ep - c_{T_1}^{4} < 0
int compute_validation_T1(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

// (A.11) in the paper.
// TODO(Javi): Fill the formula
int compute_validation_T2(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

int compute_nu(arb_ptr out, void* params, slong prec);
int compute_mu(arb_ptr out, void* params, slong prec);

#endif
