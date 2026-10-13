/* derivative_bounds.c: derivative and curvature bounds on the compact
 * interval (stage derivative_compact).
 *
 * Covers [eps, pi] (eps = 0.05) recursively with mean-value cell enclosures
 * of F_1 = R', F_2 = x R'' (R = V/(2U) - U/2), the coefficient functions
 * alpha_1, alpha_2, alpha_f and the curvature x u_0'', and compares them
 * with the targets of the summary table of Appendix B. The cover record
 * (leaves, nodes) is checked by run.py. With --reject-test one target is
 * replaced by an impossible value and the program must return UNRESOLVED;
 * --sample prints diagnostics only.
 */
#include "derivative_bounds.h"
#include <stdlib.h>
#include <time.h>
#define P CAP_PREC
static u_params_struct up;
static long nodes, leaves;
static arb_struct maxima[6];
static const char *names[] = {"F1", "F2", "alpha1", "alpha2", "alphaf", "x*u0''"};
static const char *targets[] = {"0.0000005", "0.00001", "0.32", "1.383", "0.414", "-0.089"};
int derivative_compute_dddu0(arb_ptr out, const arb_t inp, void *params, slong prec) {
  u_params_struct *up = params;
  arb_t a0kterm, a1kterm, biterm;
  arb_t aux;
  arb_t clausen;
  arb_t argument;

  int i;

  arb_init(a0kterm);
  arb_init(a1kterm);
  arb_init(biterm);
  arb_zero(a0kterm);
  arb_zero(a1kterm);
  arb_zero(biterm);

  arb_init(aux);
  arb_init(clausen);
  arb_init(argument);

  for (i = 0; i < up->N0; i++) {
    arb_mul_si(argument, up->p0, i, prec);
    arb_set_str(aux, "-1.5", prec);
    arb_add(argument, argument, aux, prec);
    arb_clausen_sin_direct(clausen, argument, inp, prec);
    arb_mul(clausen, clausen, up->a0k + i, prec);
    arb_neg(clausen, clausen);
    arb_add(a0kterm, a0kterm, clausen, prec);
  }

  for (i = 0; i < up->N1; i++) {
    arb_mul_si(argument, up->p0, i, prec);
    arb_add(argument, argument, up->p1, prec);
    arb_set_str(aux, "-1.5", prec);
    arb_add(argument, argument, aux, prec);
    arb_clausen_sin_direct(clausen, argument, inp, prec);
    arb_mul(clausen, clausen, up->a1k + i, prec);
    arb_neg(clausen, clausen);
    arb_add(a1kterm, a1kterm, clausen, prec);
  }

  for (i = 1; i <= up->N2; i++) {
    arb_mul_si(aux, inp, i, prec);
    arb_sin(aux, aux, prec);
    arb_mul_si(aux, aux, i, prec);
    arb_mul_si(aux, aux, i, prec);
    arb_mul_si(aux, aux, i, prec);
    arb_mul(aux, aux, up->bi + i - 1, prec);
    arb_add(biterm, biterm, aux, prec);
  }

  arb_add(out, a0kterm, a1kterm, prec);
  arb_add(out, out, biterm, prec);

  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);

  arb_clear(aux);
  arb_clear(clausen);
  arb_clear(argument);
  return 0;
}

int derivative_compute_dddLu0(arb_ptr out, const arb_t inp, void *params, slong prec) {
  u_params_struct *up = params;
  arb_t a0kterm, a1kterm, biterm;
  arb_t aux, aux2;
  arb_t clausen;
  arb_t zeta;
  arb_t argument;
  arb_t error_tail;
  arb_t pm_one;

  int i, m;

  arb_init(a0kterm);
  arb_init(a1kterm);
  arb_init(biterm);
  arb_zero(a0kterm);
  arb_zero(a1kterm);
  arb_zero(biterm);

  arb_init(aux);
  arb_init(aux2);
  arb_init(clausen);
  arb_init(zeta);
  arb_init(argument);
  arb_init(error_tail);
  arb_init(pm_one);
  arb_zero_pm_one(pm_one);

  // Main term
  for (i = 0; i < up->N0; i++) {
    arb_mul_si(argument, up->p0, i, prec);
    arb_set_str(aux, "-1", prec);
    arb_add(argument, argument, aux, prec);
    arb_clausen_sin_direct(clausen, argument, inp, prec);
    arb_mul(aux, clausen, up->a0k + i, prec);
    arb_neg(aux, aux);
    arb_add(a0kterm, a0kterm, aux, prec);
  }

  // Main term
  for (i = 0; i < up->N1; i++) {
    arb_mul_si(argument, up->p0, i, prec);
    arb_add(argument, argument, up->p1, prec);
    arb_set_str(aux, "-1", prec);
    arb_add(argument, argument, aux, prec);
    arb_clausen_sin_direct(clausen, argument, inp, prec);
    arb_mul(aux, clausen, up->a1k + i, prec);
    arb_neg(aux, aux);
    arb_add(a1kterm, a1kterm, aux, prec);
  }

  // Correction N0
  for (i = 0; i < up->N0; i++) {
    arb_mul_si(argument, up->p0, i, prec);
    arb_set_str(aux, "-1", prec);
    arb_add(argument, argument, aux, prec);
    arb_neg(argument, argument);

    compute_tail_error(error_tail, argument, TANH_MAX + 1, prec);
    arb_mul(error_tail, error_tail, up->a0k + i, prec);
    for (m = 1; m <= TANH_MAX; m++) {

      arb_set_si(aux2, m);
      arb_tanh(aux2, aux2, prec);
      arb_sqrt(aux2, aux2, prec);
      arb_sub_si(aux2, aux2, 1, prec);

      arb_mul_si(aux, inp, m, prec);
      arb_sin(aux, aux, prec);
      arb_neg(aux, aux);

      arb_mul(aux, aux, aux2, prec);

      arb_set_si(aux2, m);
      arb_pow(aux2, aux2, argument, prec);

      arb_mul(aux, aux, aux2, prec);
      arb_mul(aux, aux, up->a0k + i, prec);
      arb_add(a0kterm, a0kterm, aux, prec);
    }
    // Tail N0
    arb_mul(error_tail, error_tail, pm_one, prec);
    arb_add(a0kterm, a0kterm, error_tail, prec);
  }

  // Correction N1
  for (i = 0; i < up->N1; i++) {
    arb_mul_si(argument, up->p0, i, prec);
    arb_add(argument, argument, up->p1, prec);
    arb_set_str(aux, "-1", prec);
    arb_add(argument, argument, aux, prec);
    arb_neg(argument, argument);

    compute_tail_error(error_tail, argument, TANH_MAX + 1, prec);
    arb_mul(error_tail, error_tail, up->a1k + i, prec);
    for (m = 1; m <= TANH_MAX; m++) {

      arb_set_si(aux2, m);
      arb_tanh(aux2, aux2, prec);
      arb_sqrt(aux2, aux2, prec);
      arb_sub_si(aux2, aux2, 1, prec);

      arb_mul_si(aux, inp, m, prec);
      arb_sin(aux, aux, prec);
      arb_neg(aux, aux);

      arb_mul(aux, aux, aux2, prec);

      arb_set_si(aux2, m);
      arb_pow(aux2, aux2, argument, prec);

      arb_mul(aux, aux, aux2, prec);
      arb_mul(aux, aux, up->a1k + i, prec);
      arb_add(a1kterm, a1kterm, aux, prec);
    }
    // Tail N1
    arb_mul(error_tail, error_tail, pm_one, prec);
    arb_add(a1kterm, a1kterm, error_tail, prec);
  }

  for (i = 1; i <= up->N2; i++) {
    // b_k*(cos(k*x)-1)
    arb_mul_si(aux, inp, i, prec);
    arb_sin(aux, aux, prec);
    //    arb_neg(aux,aux);
    arb_mul(aux, aux, up->bi + i - 1, prec);

    //sqrt(k*tanh(k))*k
    arb_set_si(aux2, i);
    arb_tanh(aux2, aux2, prec);
    arb_mul_si(aux2, aux2, i, prec);
    arb_sqrt(aux2, aux2, prec);
    arb_mul_si(aux2, aux2, i, prec);
    arb_mul_si(aux2, aux2, i, prec);

    arb_mul(aux, aux, aux2, prec);
    arb_add(biterm, biterm, aux, prec);
  }

  arb_add(out, a0kterm, a1kterm, prec);
  arb_add(out, out, biterm, prec);

  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);

  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(clausen);
  arb_clear(zeta);
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pm_one);
  return 0;
}

void derivative_defect_jets(arb_ptr out, const arb_t x, u_params_struct *u, int order, slong p) {
  if (order < 0 || order > 3) {
    arb_indeterminate(out);
    return;
  }
  arb_ptr a = _arb_vec_init(4), b = _arb_vec_init(4), q = _arb_vec_init(4);
  arb_t t;
  arb_init(t);
  compute_u0(a, x, u, p);
  compute_Lu0(b, x, u, p);
  if (order >= 1) {
    compute_du0(a + 1, x, u, p);
    compute_dLu0(b + 1, x, u, p);
  }
  if (order >= 2) {
    compute_ddu0(a + 2, x, u, p);
    compute_ddLu0(b + 2, x, u, p);
    arb_mul_2exp_si(a + 2, a + 2, -1);
    arb_mul_2exp_si(b + 2, b + 2, -1);
  }
  if (order >= 3) {
    derivative_compute_dddu0(a + 3, x, u, p);
    derivative_compute_dddLu0(b + 3, x, u, p);
    arb_div_ui(a + 3, a + 3, 6, p);
    arb_div_ui(b + 3, b + 3, 6, p);
  }
  for (int i = 0; i <= order; i++) {
    arb_set(t, b + i);
    for (int j = 1; j <= i; j++)
      arb_submul(t, a + j, q + i - j, p);
    arb_div(q + i, t, a, p);
    arb_sub(out + i, q + i, a + i, p);
    arb_mul_2exp_si(out + i, out + i, -1);
  }
  arb_clear(t);
  _arb_vec_clear(a, 4);
  _arb_vec_clear(b, 4);
  _arb_vec_clear(q, 4);
}
static void report(const char *s, const arb_t a) {
  printf("%s = ", s);
  arb_printn(a, 23, 0);
  puts("");
  char *t = arb_dump_str(a);
  printf("ARB_DUMP %s: %s\n", s, t);
  flint_free(t);
  fflush(stdout);
}
static void coeffs(arb_ptr out, const arb_t x) {
  compute_alpha1(out, x, &up, P);
  compute_alpha2(out + 1, x, &up, P);
  compute_alphaf(out + 2, x, &up, P);
  compute_convexity(out + 3, x, &up, P);
}
static int compact(const arb_t l, const arb_t r, int depth) {
  arb_t x, m, rad, t, goal;
  arb_init(x);
  arb_init(m);
  arb_init(rad);
  arb_init(t);
  arb_init(goal);
  arb_ptr v = _arb_vec_init(6), c = _arb_vec_init(4), j = _arb_vec_init(4);
  arb_union(x, l, r, P);
  arb_add(m, l, r, P);
  arb_mul_2exp_si(m, m, -1);
  arb_sub(rad, x, m, P);
  derivative_defect_jets(c, m, &up, 2, P);
  derivative_defect_jets(j, x, &up, 3, P);
  /* F1(mid)+delta R''(X), F2(mid)+delta (R''(X)+X R'''(X)). */
  arb_mul_ui(t, j + 2, 2, P);
  arb_addmul(c + 1, rad, t, P);
  arb_abs(v, c + 1);
  arb_mul(v + 1, c + 2, m, P);
  arb_mul_ui(v + 1, v + 1, 2, P);
  arb_mul(t, j + 3, x, P);
  arb_mul_ui(t, t, 6, P);
  arb_addmul_ui(t, j + 2, 2, P);
  arb_addmul(v + 1, rad, t, P);
  arb_abs(v + 1, v + 1);
  coeffs(v + 2, x);
  int ok = 1;
  nodes++;
  for (int i = 0; i < 6; i++) {
    arb_set_str(goal, targets[i], P);
    if (!arb_is_finite(v + i) || !arb_lt(v + i, goal))
      ok = 0;
  }
  if (ok) {
    leaves++;
    for (int i = 0; i < 6; i++)
      arb_max(maxima + i, maxima + i, v + i, P);
    if (leaves % 2000 == 0) {
      printf("leaves=%ld nodes=%ld x=", leaves, nodes);
      arb_printn(m, 8, 0);
      puts("");
      fflush(stdout);
    }
  } else if (depth < 32 && nodes < 2000000) {
    ok = compact(l, m, depth + 1);
    if (ok)
      ok = compact(m, r, depth + 1);
  } else {
    report("UNRESOLVED cell", x);
    for (int i = 0; i < 6; i++)
      report(names[i], v + i);
  }
  _arb_vec_clear(v, 6);
  _arb_vec_clear(c, 4);
  _arb_vec_clear(j, 4);
  arb_clear(x);
  arb_clear(m);
  arb_clear(rad);
  arb_clear(t);
  arb_clear(goal);
  return ok;
}
int main(int argc, char **argv) {
  populate_parameters(&up, P);
  arb_t x, pi;
  arb_init(x);
  arb_init(pi);
  arb_const_pi(pi, P);
  arb_ptr v = _arb_vec_init(6), j = _arb_vec_init(4);
  for (int i = 0; i < 6; i++) {
    arb_init(maxima + i);
    arb_set_si(maxima + i, -100);
  }
  if (argc > 1 && !strcmp(argv[1], "--sample")) {
    for (int k = 1; k <= 256; k++) {
      arb_mul_ui(x, pi, k, P);
      arb_div_ui(x, x, 256, P);
      derivative_defect_jets(j, x, &up, 2, P);
      arb_set(v, j + 1);
      arb_mul(v + 1, j + 2, x, P);
      arb_mul_ui(v + 1, v + 1, 2, P);
      coeffs(v + 2, x);
      report("x", x);
      for (int i = 0; i < 6; i++)
        report(names[i], v + i);
    }
    return 0;
  }
  if (argc > 2 && !strcmp(argv[2], "--reject-test"))
    targets[0] = "0";
  arb_set_str(x, argc > 1 ? argv[1] : "0.05", P);
  int ok = compact(x, pi, 0);
  for (int i = 0; i < 6; i++)
    report(names[i], maxima + i);
  printf("DERIVATIVE COMPACT %s leaves=%ld nodes=%ld\n", ok ? "PASS" : "UNRESOLVED", leaves, nodes);
  return ok ? 0 : 2;
}
