/* cusp.c: the cusp expansions shared by the convexity programs; see cusp.h.
 *
 * The coefficient tables are built once from the exact approximation. The
 * Clausen tails after CUSP_M - 1 even terms and the Fourier tails of
 * Appendix B give the remainders; derivative_profile_jets() encloses
 * x^j U^(j)(x)/sqrt(x), j = 0, 1, 2, on a whole interval 0 < X <= 0.1.
 */
#include "cusp.h"
#include <stdlib.h>
#define P CAP_PREC
u_params_struct cusp_up;
cusp_term cusp_us[24], cusp_ss[512];
int cusp_nu, cusp_ns;
void cusp_exponent(arb_t a, int h, int k, int j) {
  arb_set_si(a, h);
  arb_mul_2exp_si(a, a, -1);
  arb_addmul_ui(a, cusp_up.p0, k, P);
  arb_addmul_ui(a, cusp_up.p1, j, P);
}
void cusp_falling(arb_t a, const arb_t e, int d) {
  arb_t t;
  arb_init(t);
  arb_one(a);
  for (int i = 0; i < d; i++) {
    arb_sub_ui(t, e, i, P);
    arb_mul(a, a, t, P);
  }
  arb_clear(t);
}
/* sum(1-sqrt(tanh n))*n^e, including its infinite tail. */
void cusp_smooth(arb_t out, const arb_t e) {
  arb_t t, a, b;
  arb_init(t);
  arb_init(a);
  arb_init(b);
  arb_zero(out);
  for (int n = 1; n <= CUSP_N; n++) {
    arb_set_ui(a, n);
    arb_tanh(t, a, P);
    arb_sqrt(t, t, P);
    arb_neg(t, t);
    arb_add_ui(t, t, 1, P);
    arb_pow(b, a, e, P);
    arb_addmul(out, t, b, P);
  }
  compute_tail_error(t, e, CUSP_N + 1, P);
  arb_add_error(out, t);
  arb_clear(t);
  arb_clear(a);
  arb_clear(b);
}
static void add(int h, int k, int j, const arb_t c) {
  int i;
  for (i = 0; i < cusp_ns; i++)
    if (cusp_ss[i].h == h && cusp_ss[i].k == k && cusp_ss[i].j == j)
      break;
  if (i == cusp_ns) {
    if (cusp_ns >= 512)
      abort();
    cusp_ns++;
    cusp_ss[i].h = h;
    cusp_ss[i].k = k;
    cusp_ss[i].j = j;
    arb_init(cusp_ss[i].c);
  }
  arb_add(cusp_ss[i].c, cusp_ss[i].c, c, P);
}
/* U or V regular x^(2m) coefficient, signed, including complete cusp_smooth sums. */
static void coefficient(arb_t out, int m, int v) {
  arb_t s, e, t, a, f;
  arb_init(s);
  arb_init(e);
  arb_init(t);
  arb_init(a);
  arb_init(f);
  arb_zero(out);
  arb_fac_ui(f, 2 * m, P);
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? cusp_up.N1 : cusp_up.N0); k++) {
      cusp_exponent(s, v ? 4 : 3, k, j);
      arb_sub_ui(e, s, 2 * m, P);
      arb_zeta(t, e, P);
      if (v) {
        arb_neg(e, e);
        cusp_smooth(a, e);
        arb_sub(t, t, a, P);
      }
      arb_mul(t, t, (j ? cusp_up.a1k : cusp_up.a0k) + k, P);
      if (!(m % 2))
        arb_neg(t, t);
      arb_add(out, out, t, P);
    }
  for (int n = 1; n <= cusp_up.N2; n++) {
    arb_set_ui(a, n);
    arb_pow_ui(t, a, 2 * m, P);
    if (v) {
      arb_tanh(a, a, P);
      arb_div_ui(a, a, n, P);
      arb_sqrt(a, a, P);
      arb_mul(t, t, a, P);
    }
    arb_mul(t, t, cusp_up.bi + n - 1, P);
    if (m % 2)
      arb_neg(t, t);
    arb_add(out, out, t, P);
  }
  arb_div(out, out, f, P);
  arb_clear(s);
  arb_clear(e);
  arb_clear(t);
  arb_clear(a);
  arb_clear(f);
}
/* Clausen analytic tail after m=CUSP_M-1, bounds x^d R^(d), d<=2.
 |zeta(s-k)|/k! <= 2 zeta(k+1-s)(2pi)^(s-1-k)(k+2)^2.
 The next even-cusp_term ratio is <=(x/2pi)^2((2M+4)/(2M+2))^(2+d). */
void cusp_cltail(arb_t out, const arb_t s, const arb_t x, int d) {
  arb_t a, b, c, r;
  arb_init(a);
  arb_init(b);
  arb_init(c);
  arb_init(r);
  arb_set_si(a, -2);
  arb_set_ui(b, 7);
  if (!arb_ge(s, a) || !arb_le(s, b) || !arb_is_positive(x))
    abort();
  arb_const_pi(a, P);
  arb_mul_ui(a, a, 2, P);
  arb_div(r, x, a, P);
  arb_pow_ui(out, r, 2 * CUSP_M, P);
  arb_sqr(r, r, P);
  arb_set_ui(b, 2 * CUSP_M + 4);
  arb_div_ui(b, b, 2 * CUSP_M + 2, P);
  arb_pow_ui(b, b, 2 + d, P);
  arb_mul(r, r, b, P);
  arb_neg(r, r);
  arb_add_ui(r, r, 1, P);
  if (!arb_is_positive(r))
    abort();
  arb_div(out, out, r, P);
  arb_sub_ui(b, s, 1, P);
  arb_pow(a, a, b, P);
  arb_mul(out, out, a, P);
  arb_set_ui(a, 2 * CUSP_M + 2);
  arb_pow_ui(a, a, 2 + d, P);
  arb_mul(out, out, a, P);
  arb_neg(b, s);
  arb_add_ui(b, b, 2 * CUSP_M + 1, P);
  arb_zeta(b, b, P);
  arb_mul(out, out, b, P);
  arb_mul_ui(out, out, 2, P);
  arb_clear(a);
  arb_clear(b);
  arb_clear(c);
  arb_clear(r);
}
void cusp_remainder(arb_t out, const arb_t x, int d, int v) {
  arb_t s, e, t, a, f;
  arb_init(s);
  arb_init(e);
  arb_init(t);
  arb_init(a);
  arb_init(f);
  arb_zero(out);
  arb_fac_ui(f, 2 * CUSP_M - d, P);
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? cusp_up.N1 : cusp_up.N0); k++) {
      cusp_exponent(s, v ? 4 : 3, k, j);
      cusp_cltail(t, s, x, d);
      if (v) {
        arb_neg(e, s);
        arb_add_ui(e, e, 2 * CUSP_M, P);
        cusp_smooth(a, e);
        arb_abs(a, a);
        arb_div(a, a, f, P);
        arb_pow_ui(e, x, 2 * CUSP_M, P);
        arb_addmul(t, a, e, P);
      }
      arb_abs(a, (j ? cusp_up.a1k : cusp_up.a0k) + k);
      arb_addmul(out, a, t, P);
    }
  for (int n = 1; n <= cusp_up.N2; n++) {
    arb_set_ui(a, n);
    arb_pow_ui(t, a, 2 * CUSP_M, P);
    if (v) {
      arb_tanh(a, a, P);
      arb_div_ui(a, a, n, P);
      arb_sqrt(a, a, P);
      arb_mul(t, t, a, P);
    }
    arb_abs(a, cusp_up.bi + n - 1);
    arb_mul(t, t, a, P);
    arb_div(t, t, f, P);
    arb_pow_ui(a, x, 2 * CUSP_M, P);
    arb_addmul(out, t, a, P);
  }
  arb_clear(s);
  arb_clear(e);
  arb_clear(t);
  arb_clear(a);
  arb_clear(f);
}
static int initialized;
static void build_tables(void) {
  if (initialized)
    return;
  initialized = 1;
  arb_t a;
  arb_init(a);
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? cusp_up.N1 : cusp_up.N0); k++) {
      cusp_us[cusp_nu].h = 1;
      cusp_us[cusp_nu].k = k;
      cusp_us[cusp_nu].j = j;
      arb_init(cusp_us[cusp_nu].c);
      compute_a_jk_0(cusp_us[cusp_nu].c, j, k, &cusp_up, P);
      compute_A_jk_0(a, j, k, &cusp_up, P);
      add(2, k, j, a);
      cusp_nu++;
    }
  for (int m = 1; m < CUSP_M; m++) {
    cusp_us[cusp_nu].h = 4 * m;
    cusp_us[cusp_nu].k = cusp_us[cusp_nu].j = 0;
    arb_init(cusp_us[cusp_nu].c);
    coefficient(cusp_us[cusp_nu].c, m, 0);
    coefficient(a, m, 1);
    add(4 * m, 0, 0, a);
    cusp_nu++;
  }
  for (int i = 0; i < cusp_nu; i++)
    for (int j = 0; j < cusp_nu; j++) {
      arb_mul(a, cusp_us[i].c, cusp_us[j].c, P);
      arb_neg(a, a);
      add(cusp_us[i].h + cusp_us[j].h, cusp_us[i].k + cusp_us[j].k, cusp_us[i].j + cusp_us[j].j, a);
    }
  arb_clear(a);
}
/* Fixes the approximation and builds the coefficient tables; later calls with
   the same approximation are no-ops. */
void cusp_init(const u_params_struct *params) {
  if (initialized)
    return;
  cusp_up = *params;
  build_tables();
}
/* Cusp remainder bounds of the scaled defect for d = 0, 1, 2 on (0, eps]. */
void derivative_cusp_errors(arb_ptr out, const arb_t eps, u_params_struct *params, slong prec) {
  if (prec != P || !arb_is_positive(eps))
    abort();
  cusp_init(params);
  arb_t a, e, t, pow, den;
  arb_init(a);
  arb_init(e);
  arb_init(t);
  arb_init(pow);
  arb_init(den);
  arb_set_str(t, "0.5", P);
  arb_add(t, t, cusp_up.p0, P);
  arb_pow(den, eps, t, P);
  for (int d = 0; d < 3; d++) {
    arb_zero(out + d);
    for (int i = 1; i < cusp_nu; i++) {
      cusp_exponent(e, cusp_us[i].h, cusp_us[i].k, cusp_us[i].j);
      cusp_falling(a, e, d);
      arb_mul(a, a, cusp_us[i].c, P);
      arb_abs(a, a);
      cusp_exponent(e, cusp_us[i].h - 1, cusp_us[i].k, cusp_us[i].j);
      arb_sub(e, e, cusp_up.p0, P);
      /* Exact lowest cusp_exponent p0-p0 is zero despite interval subtraction. */
      if (cusp_us[i].h == 1 && cusp_us[i].k == 1 && cusp_us[i].j == 0)
        arb_one(pow);
      else {
        if (!arb_is_positive(e))
          abort();
        arb_pow(pow, eps, e, P);
      }
      arb_addmul(out + d, a, pow, P);
    }
    cusp_remainder(t, eps, d, 0);
    arb_div(t, t, den, P);
    arb_add(out + d, out + d, t, P);
  }
  arb_clear(a);
  arb_clear(e);
  arb_clear(t);
  arb_clear(pow);
  arb_clear(den);
}
/* out[j] encloses x^j U^(j)(x)/sqrt(x), j=0,1,2, for all x in X.
 * The source's exact input is fixed at first call, as in derivative_cusp_errors.
 * Returns 1 only for a finite enclosure on 0<X<=0.1 at the working precision CAP_PREC.
 * No mutable global manuscript or existence code is used. */
int derivative_profile_jets(arb_ptr out, const arb_t X, u_params_struct *params, slong prec) {
  arb_t limit, eps, e, a, power, tail, sq;
  arb_init(limit);
  arb_init(eps);
  arb_init(e);
  arb_init(a);
  arb_init(power);
  arb_init(tail);
  arb_init(sq);
  arb_set_str(limit, "0.1", P);
  int ok = prec == P && arb_is_positive(X) && arb_le(X, limit);
  if (!ok)
    goto done;
  cusp_init(params);
  arf_t upper;
  arf_init(upper);
  arb_get_ubound_arf(upper, X, P);
  arb_set_arf(eps, upper);
  arf_clear(upper);
  arb_sqrt(sq, eps, P);
  for (int d = 0; d < 3; ++d) {
    arb_zero(out + d);
    for (int i = 0; i < cusp_nu; ++i) {
      cusp_exponent(e, cusp_us[i].h, cusp_us[i].k, cusp_us[i].j);
      cusp_falling(a, e, d);
      arb_mul(a, a, cusp_us[i].c, P);
      cusp_exponent(e, cusp_us[i].h - 1, cusp_us[i].k, cusp_us[i].j);
      arb_pow(power, X, e, P);
      arb_addmul(out + d, a, power, P);
    }
    cusp_remainder(tail, eps, d, 0);
    /* Tail/sqrt(x) is a positive increasing majorant. Its endpoint
     * evaluation bounds the whole X cell without division by sqrt(X). */
    arb_div(tail, tail, sq, P);
    arb_add_error(out + d, tail);
    ok = ok && arb_is_finite(out + d);
  }
done:
  if (!ok)
    for (int d = 0; d < 3; ++d)
      arb_indeterminate(out + d);
  arb_clear(limit);
  arb_clear(eps);
  arb_clear(e);
  arb_clear(a);
  arb_clear(power);
  arb_clear(tail);
  arb_clear(sq);
  return ok;
}
