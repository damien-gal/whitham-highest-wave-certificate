/* derivative_cusp.c: derivative and curvature bounds on the cusp interval
 * (stage derivative_cusp) and the cusp jets used by the other convexity
 * programs.
 *
 * On 0 < x <= eps (eps = 0.05) the scaled quantities U = u_0/sqrt(x) and the
 * residual of the equation are finite sums of powers x^(h/2 + k p_0 + j p_1)
 * with explicit remainders (Clausen tails after m = M-1 even terms and the
 * Fourier tails of Appendix B). The program bounds |F_1|, |F_2|, the
 * coefficient functions alpha_1, alpha_2, alpha_f and -x u_0'' against the
 * targets of the summary table of Appendix B, together with the printed
 * constants C < 0.035061, |x f/(2 u_0)| < 0.114 and x u_0'' < -0.70.
 * Compiled with DERIVATIVE_CUSP_LIBRARY, the file provides only
 * derivative_cusp_errors(), the remainder bounds used by derivative_profile.c.
 */
#include "constant_functions.h"
#include <stdlib.h>
#define P CAP_PREC
#define M 7
#define N 128
static u_params_struct up;
typedef struct {
  int h, k, j;
  arb_t c;
} term;
static term us[24], ss[512];
static int nu, ns;
#ifndef DERIVATIVE_CUSP_LIBRARY
static void show(const char *n, const arb_t x) {
  printf("%s = ", n);
  arb_printn(x, 25, 0);
  puts("");
  char *s = arb_dump_str(x);
  printf("ARB_DUMP %s: %s\n", n, s);
  flint_free(s);
}
#endif
static void exponent(arb_t a, int h, int k, int j) {
  arb_set_si(a, h);
  arb_mul_2exp_si(a, a, -1);
  arb_addmul_ui(a, up.p0, k, P);
  arb_addmul_ui(a, up.p1, j, P);
}
static void falling(arb_t a, const arb_t e, int d) {
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
static void smooth(arb_t out, const arb_t e) {
  arb_t t, a, b;
  arb_init(t);
  arb_init(a);
  arb_init(b);
  arb_zero(out);
  for (int n = 1; n <= N; n++) {
    arb_set_ui(a, n);
    arb_tanh(t, a, P);
    arb_sqrt(t, t, P);
    arb_neg(t, t);
    arb_add_ui(t, t, 1, P);
    arb_pow(b, a, e, P);
    arb_addmul(out, t, b, P);
  }
  compute_tail_error(t, e, N + 1, P);
  arb_add_error(out, t);
  arb_clear(t);
  arb_clear(a);
  arb_clear(b);
}
static void add(int h, int k, int j, const arb_t c) {
  int i;
  for (i = 0; i < ns; i++)
    if (ss[i].h == h && ss[i].k == k && ss[i].j == j)
      break;
  if (i == ns) {
    if (ns >= 512)
      abort();
    ns++;
    ss[i].h = h;
    ss[i].k = k;
    ss[i].j = j;
    arb_init(ss[i].c);
  }
  arb_add(ss[i].c, ss[i].c, c, P);
}
/* U or V regular x^(2m) coefficient, signed, including complete smooth sums. */
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
    for (int k = 0; k < (j ? up.N1 : up.N0); k++) {
      exponent(s, v ? 4 : 3, k, j);
      arb_sub_ui(e, s, 2 * m, P);
      arb_zeta(t, e, P);
      if (v) {
        arb_neg(e, e);
        smooth(a, e);
        arb_sub(t, t, a, P);
      }
      arb_mul(t, t, (j ? up.a1k : up.a0k) + k, P);
      if (!(m % 2))
        arb_neg(t, t);
      arb_add(out, out, t, P);
    }
  for (int n = 1; n <= up.N2; n++) {
    arb_set_ui(a, n);
    arb_pow_ui(t, a, 2 * m, P);
    if (v) {
      arb_tanh(a, a, P);
      arb_div_ui(a, a, n, P);
      arb_sqrt(a, a, P);
      arb_mul(t, t, a, P);
    }
    arb_mul(t, t, up.bi + n - 1, P);
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
/* Clausen analytic tail after m=M-1, bounds x^d R^(d), d<=2.
 |zeta(s-k)|/k! <= 2 zeta(k+1-s)(2pi)^(s-1-k)(k+2)^2.
 The next even-term ratio is <=(x/2pi)^2((2M+4)/(2M+2))^(2+d). */
static void cltail(arb_t out, const arb_t s, const arb_t x, int d) {
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
  arb_pow_ui(out, r, 2 * M, P);
  arb_sqr(r, r, P);
  arb_set_ui(b, 2 * M + 4);
  arb_div_ui(b, b, 2 * M + 2, P);
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
  arb_set_ui(a, 2 * M + 2);
  arb_pow_ui(a, a, 2 + d, P);
  arb_mul(out, out, a, P);
  arb_neg(b, s);
  arb_add_ui(b, b, 2 * M + 1, P);
  arb_zeta(b, b, P);
  arb_mul(out, out, b, P);
  arb_mul_ui(out, out, 2, P);
  arb_clear(a);
  arb_clear(b);
  arb_clear(c);
  arb_clear(r);
}
static void cusp_remainder(arb_t out, const arb_t x, int d, int v) {
  arb_t s, e, t, a, f;
  arb_init(s);
  arb_init(e);
  arb_init(t);
  arb_init(a);
  arb_init(f);
  arb_zero(out);
  arb_fac_ui(f, 2 * M - d, P);
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? up.N1 : up.N0); k++) {
      exponent(s, v ? 4 : 3, k, j);
      cltail(t, s, x, d);
      if (v) {
        arb_neg(e, s);
        arb_add_ui(e, e, 2 * M, P);
        smooth(a, e);
        arb_abs(a, a);
        arb_div(a, a, f, P);
        arb_pow_ui(e, x, 2 * M, P);
        arb_addmul(t, a, e, P);
      }
      arb_abs(a, (j ? up.a1k : up.a0k) + k);
      arb_addmul(out, a, t, P);
    }
  for (int n = 1; n <= up.N2; n++) {
    arb_set_ui(a, n);
    arb_pow_ui(t, a, 2 * M, P);
    if (v) {
      arb_tanh(a, a, P);
      arb_div_ui(a, a, n, P);
      arb_sqrt(a, a, P);
      arb_mul(t, t, a, P);
    }
    arb_abs(a, up.bi + n - 1);
    arb_mul(t, t, a, P);
    arb_div(t, t, f, P);
    arb_pow_ui(a, x, 2 * M, P);
    arb_addmul(out, t, a, P);
  }
  arb_clear(s);
  arb_clear(e);
  arb_clear(t);
  arb_clear(a);
  arb_clear(f);
}
#ifndef DERIVATIVE_CUSP_LIBRARY
/* Signed interval for x^d P^(d)/x^(shift/2) throughout (0,eps]. */
static void polynomial(arb_t out, term *ts, int n, int shift, int d, const arb_t eps, int skip) {
  arb_t a, e, t, b, z;
  arb_init(a);
  arb_init(e);
  arb_init(t);
  arb_init(b);
  arb_init(z);
  arb_zero(out);
  for (int i = 0; i < n; i++) {
    if (skip && ts[i].h == 2 && ts[i].k == 0 && ts[i].j == 0)
      continue;
    exponent(e, ts[i].h, ts[i].k, ts[i].j);
    falling(a, e, d);
    arb_mul(a, a, ts[i].c, P);
    exponent(e, ts[i].h - shift, ts[i].k, ts[i].j);
    if (arb_is_zero(e))
      arb_one(t);
    else {
      if (!arb_is_positive(e))
        abort();
      arb_pow(b, eps, e, P);
      arb_union(t, z, b, P);
    }
    arb_addmul(out, a, t, P);
  }
  arb_clear(a);
  arb_clear(e);
  arb_clear(t);
  arb_clear(b);
  arb_clear(z);
}
#endif
static int initialized;
static void build_tables(void) {
  if (initialized)
    return;
  initialized = 1;
  arb_t a;
  arb_init(a);
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? up.N1 : up.N0); k++) {
      us[nu].h = 1;
      us[nu].k = k;
      us[nu].j = j;
      arb_init(us[nu].c);
      compute_a_jk_0(us[nu].c, j, k, &up, P);
      compute_A_jk_0(a, j, k, &up, P);
      add(2, k, j, a);
      nu++;
    }
  for (int m = 1; m < M; m++) {
    us[nu].h = 4 * m;
    us[nu].k = us[nu].j = 0;
    arb_init(us[nu].c);
    coefficient(us[nu].c, m, 0);
    coefficient(a, m, 1);
    add(4 * m, 0, 0, a);
    nu++;
  }
  for (int i = 0; i < nu; i++)
    for (int j = 0; j < nu; j++) {
      arb_mul(a, us[i].c, us[j].c, P);
      arb_neg(a, a);
      add(us[i].h + us[j].h, us[i].k + us[j].k, us[i].j + us[j].j, a);
    }
  arb_clear(a);
}
void derivative_cusp_errors(arb_ptr out, const arb_t eps, u_params_struct *params, slong prec) {
  if (prec != P || !arb_is_positive(eps))
    abort();
  if (!initialized) {
    up = *params;
    build_tables();
  }
  arb_t a, e, t, pow, den;
  arb_init(a);
  arb_init(e);
  arb_init(t);
  arb_init(pow);
  arb_init(den);
  arb_set_str(t, "0.5", P);
  arb_add(t, t, up.p0, P);
  arb_pow(den, eps, t, P);
  for (int d = 0; d < 3; d++) {
    arb_zero(out + d);
    for (int i = 1; i < nu; i++) {
      exponent(e, us[i].h, us[i].k, us[i].j);
      falling(a, e, d);
      arb_mul(a, a, us[i].c, P);
      arb_abs(a, a);
      exponent(e, us[i].h - 1, us[i].k, us[i].j);
      arb_sub(e, e, up.p0, P);
      /* Exact lowest exponent p0-p0 is zero despite interval subtraction. */
      if (us[i].h == 1 && us[i].k == 1 && us[i].j == 0)
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
#ifndef DERIVATIVE_CUSP_LIBRARY
int main(int argc, char **argv) {
  populate_parameters(&up, P);
  arb_t eps, t, a, b, c, e, lambda, sqrtx;
  arb_init(eps);
  arb_init(t);
  arb_init(a);
  arb_init(b);
  arb_init(c);
  arb_init(e);
  arb_init(lambda);
  arb_init(sqrtx);
  arb_set_str(eps, argc > 1 ? argv[1] : "0.05", P);
  arb_sqrt(sqrtx, eps, P);
  arb_set_str(a, "0.25", P);
  if (!arb_equal(a, up.a0k))
    return 2;
  build_tables();
  arb_ptr U = _arb_vec_init(3), B = _arb_vec_init(3), Eu = _arb_vec_init(3), Ev = _arb_vec_init(3),
          Pu = _arb_vec_init(3), vals = _arb_vec_init(7);
  for (int d = 0; d < 3; d++) {
    polynomial(Pu + d, us, nu, 1, d, eps, 0);
    arb_set(U + d, Pu + d);
    cusp_remainder(Eu + d, eps, d, 0);
    arb_div(t, Eu + d, sqrtx, P);
    arb_add_error(U + d, t);
    cusp_remainder(Ev + d, eps, d, 1);
    polynomial(B + d, ss, ns, 3, d, eps, 1);
    arb_mul(t, eps, sqrtx, P);
    arb_div(e, Ev + d, t, P);
    for (int j = 0; j <= d; j++) {
      int choose = d == 2 && j == 1 ? 2 : 1;
      arb_abs(a, Pu + j);
      arb_mul(a, a, Eu + d - j, P);
      arb_mul_ui(a, a, 2 * choose, P);
      arb_div(a, a, eps, P);
      arb_add(e, e, a, P);
      arb_mul(a, Eu + j, Eu + d - j, P);
      arb_mul_ui(a, a, choose, P);
      arb_div(a, a, t, P);
      arb_add(e, e, a, P);
    }
    arb_add_error(B + d, e);
  }
  /* R=S/(2U), F1=R', F2=xR''. */
  arb_div(a, U + 1, U, P);
  arb_mul(t, B, a, P);
  arb_sub(vals, B + 1, t, P);
  arb_div(vals, vals, U, P);
  arb_mul_2exp_si(vals, vals, -1);
  arb_abs(vals, vals);
  arb_mul(t, B + 1, a, P);
  arb_mul_ui(t, t, 2, P);
  arb_sub(vals + 1, B + 2, t, P);
  arb_div(b, U + 2, U, P);
  arb_submul(vals + 1, B, b, P);
  arb_sqr(c, a, P);
  arb_mul_ui(c, c, 2, P);
  arb_addmul(vals + 1, B, c, P);
  arb_div(vals + 1, vals + 1, U, P);
  arb_mul_2exp_si(vals + 1, vals + 1, -1);
  arb_abs(vals + 1, vals + 1);
  arb_sqr(t, U, P);
  arb_mul_ui(t, t, 2, P);
  arb_mul(vals + 2, U + 1, sqrtx, P);
  arb_div(vals + 2, vals + 2, t, P);
  arb_abs(vals + 2, vals + 2);
  arb_sqr(c, U + 1, P);
  arb_mul_ui(c, c, 2, P);
  arb_div(c, c, U, P);
  arb_sub(c, U + 2, c, P);
  arb_mul(vals + 3, c, sqrtx, P);
  arb_div(vals + 3, vals + 3, t, P);
  arb_abs(vals + 3, vals + 3);
  arb_div(vals + 5, U + 2, sqrtx, P); /* Since negative numerator, valid uniform upper bound. */
  /* Full-mode alpha_f: K=kappa x^-1/2+g0+O(C x²). */
  arb_t pi, kappa, g0, C;
  arb_init(pi);
  arb_init(kappa);
  arb_init(g0);
  arb_init(C);
  arb_const_pi(pi, P);
  arb_mul_ui(kappa, pi, 2, P);
  arb_sqrt(kappa, kappa, P);
  arb_inv(kappa, kappa, P);
  arb_set_str(e, "0.5", P);
  arb_zeta(g0, e, P);
  arb_set_str(e, "-0.5", P);
  smooth(t, e);
  arb_sub(g0, g0, t, P);
  arb_set_str(t, "0.5", P);
  arb_add(g0, g0, t, P);
  arb_div(g0, g0, pi, P);
  arb_zero(C);
  for (int m = 1; m < M; m++) {
    arb_set_str(e, "0.5", P);
    arb_sub_ui(e, e, 2 * m, P);
    arb_zeta(t, e, P);
    arb_abs(t, t);
    arb_fac_ui(e, 2 * m, P);
    arb_div(t, t, e, P);
    arb_pow_ui(e, eps, 2 * m - 2, P);
    arb_addmul(C, t, e, P);
  }
  arb_set_str(e, "0.5", P);
  cltail(t, e, eps, 0);
  arb_sqr(e, eps, P);
  arb_div(t, t, e, P);
  arb_add(C, C, t, P);
  arb_set_str(e, "1.5", P);
  smooth(t, e);
  arb_abs(t, t);
  arb_mul_2exp_si(t, t, -1);
  arb_add(C, C, t, P);
  arb_div(C, C, pi, P);
  /* leading factor = -(4/3)t+(8/3)t²-(4/3)h, t=a-1/2,h=b+1/4 */
  arb_set_str(t, "0.5", P);
  arb_sub(t, a, t, P);
  arb_set_str(e, "0.25", P);
  arb_add(e, b, e, P);
  arb_sqr(c, t, P);
  arb_mul_ui(c, c, 8, P);
  arb_submul_ui(c, t, 4, P);
  arb_submul_ui(c, e, 4, P);
  arb_div_ui(c, c, 3, P);
  arb_mul(c, c, kappa, P);
  arb_abs(vals + 4, c);
  arb_sqr(c, a, P);
  arb_mul_ui(c, c, 2, P);
  arb_sub(c, c, b, P); /* c=2a²-b */
  arb_mul_2exp_si(t, c, -1);
  arb_submul_ui(t, a, 2, P);
  arb_add_ui(t, t, 1, P);
  arb_mul(t, t, g0, P);
  arb_abs(t, t);
  arb_mul(t, t, sqrtx, P);
  arb_add(vals + 4, vals + 4, t, P);
  arb_abs(t, a);
  arb_mul_ui(t, t, 2, P);
  arb_div_ui(t, t, 3, P);
  arb_add_ui(t, t, 1, P);
  arb_abs(e, c);
  arb_div_ui(e, e, 12, P);
  arb_add(t, t, e, P);
  arb_mul(t, t, C, P);
  arb_sqr(e, eps, P);
  arb_mul(t, t, e, P);
  arb_mul(t, t, sqrtx, P);
  arb_add(vals + 4, vals + 4, t, P);
  arb_div(vals + 4, vals + 4, U, P);
  const char *names[] = {"cusp F1",     "cusp F2",     "cusp alpha1",
                         "cusp alpha2", "cusp alphaf", "cusp x*u0''"};
  const char *targets[] = {"0.0000005", "0.00001", "0.32", "1.383", "0.414", "-0.089"};
  int ok = arb_is_positive(U);
  for (int d = 0; d < 3; d++) {
    char s[40];
    sprintf(s, "scaled U derivative %d", d);
    show(s, U + d);
    sprintf(s, "scaled residual derivative %d", d);
    show(s, B + d);
  }
  show("kernel regular g0", g0);
  show("kernel regular O(x²) bound", C);
  for (int i = 0; i < 6; i++) {
    show(names[i], vals + i);
    arb_set_str(t, targets[i], P);
    if (!arb_is_finite(vals + i) || !arb_lt(vals + i, t))
      ok = 0;
  }
  /* Printed in Appendix B: C < 0.035061, |xf/(2u0)| < 0.114, x*u0'' < -0.70 on the cusp. */
  arb_set_str(t, "0.035061", P);
  if (!arb_lt(C, t))
    ok = 0;
  arb_set_str(t, "0.114", P);
  if (!arb_lt(vals + 4, t))
    ok = 0;
  arb_set_str(t, "-0.70", P);
  if (!arb_lt(vals + 5, t))
    ok = 0;
  printf("DERIVATIVE CUSP %s\n", ok ? "PASS" : "UNRESOLVED");
  return ok ? 0 : 2;
}

#endif
