/* constant_functions.h: coefficients, remainders and scalar constants of the
 * cusp expansions of u_0 and L u_0 (Section 5 and Appendix B of the paper).
 *
 * Near x = 0 the Clausen series of u_0 and of L u_0 split into singular
 * powers x^(1/2 + k p_0 + j p_1) and even powers x^(2m). The routines below
 * enclose their coefficients, the remainders after M even terms and the
 * scalar constants of the lemma on u_0 near the cusp, of the local T_0
 * criterion and of the wave speed. Every routine stores a ball in out; the
 * integer return value is not used by any program (see basic_functions.h).
 */
#ifndef CONSTANT_FUNCTIONS_H_INCLUDED
#define CONSTANT_FUNCTIONS_H_INCLUDED

#include "basic_functions.h"

/* Coefficient of the singular power x^(1/2 + k p_0 + j p_1) (j = 0, 1) in the
   expansion of u_0 (a) and of L u_0 (A). Returns 0 for an index outside the
   expansion, 1 otherwise. */
int compute_a_jk_0(arb_ptr out, int j, int k, void* params, slong prec);
int compute_A_jk_0(arb_ptr out, int j, int k, void* params, slong prec);

/* Coefficient of the even power x^(2m) in the expansion of u_0 (a) and of
   L u_0 (A), including the complete sums over the tanh weights. */
int compute_a_m_0(arb_ptr out, int m, void* params, slong prec);
int compute_A_m_0(arb_ptr out, int m, void* params, slong prec);

/* The coefficients of x^2 and x^4 of u_0 (m = 1 and m = 2 above). */
int compute_b_0_0(arb_ptr out, void* params, slong prec);
int compute_b_1_0(arb_ptr out, void* params, slong prec);

/* Remainder bounds for the expansions of u_0 and L u_0 truncated after M even
   terms, valid for 0 < x <= inp (Appendix B). */
int compute_E_x_u0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);
int compute_E_x_Lu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec);

/* Constants of the lemma on u_0 near the cusp, for 0 < x <= eps with M terms:
   c_hat_plus/minus are the upper and lower hat constants of the normalized
   profile, c_u0_hat_plus/minus the constants c_hat/(lambda -+ c_hat eps^p0)
   (printed as 0.24 for the upper one), and c_p0 = c_hat_plus is the constant
   printed as 0.142. */
int compute_c_hat_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_hat_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_u0_hat_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_u0_hat_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec);
int compute_c_p0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

/* The constants c_{T_0}^1(eps) and c_{T_0}^2(eps) of the local T_0 estimate
   (Section 6), and the left-hand side minus the right-hand side of the
   criterion (ineqT0) with C_B replaced by 1 and c_{p_0} by c_p0(eps, M). */
int compute_cT0_1(arb_ptr out, const arb_t eps, void* params, slong prec);
int compute_cT0_2(arb_ptr out, const arb_t eps, void* params, slong prec);
int compute_validation_T0(arb_ptr out, const arb_t eps, int M, void* params, slong prec);

/* nu = 4 L u_0(0) = 8 int_0^pi K u_0 dy, with the symmetric error 4 pi eps0
   of the wave-speed proposition. */
int compute_nu(arb_ptr out, const arb_t eps0, void* params, slong prec);

#endif
