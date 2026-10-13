/* verify.c: the bound on the defect F_0 (stage existence).
 *
 * Usage: verify <eps> compact [--reject-test]. On (0, eps] (eps = 0.05) the
 * defect is bounded from the cusp expansions of u_0 and L u_0 with the
 * remainder bounds of constant_functions.c; on [eps, pi] a recursive cover
 * with mean-value cell enclosures certifies |F_0| <= 1e-7 and prints the
 * cover record checked by run.py. The constants printed in Appendix B
 * (lambda_- > 0.604 and the cusp defect bound 1.55e-8) are asserted. With
 * --reject-test an impossible target is used and the program must return
 * UNRESOLVED.
 */
#include "constant_functions.h"
#include "complex_functions.h"
#include <stdlib.h>
#include <time.h>
#define PREC CAP_PREC
static slong p = PREC;
static u_params_struct up;
static void show(const char *s, const arb_t x) {
  printf("%s = ", s);
  arb_printn(x, 25, ARB_STR_MORE);
  puts("");
  char *dump = arb_dump_str(x);
  printf("ARB_DUMP %s: %s\n", s, dump);
  flint_free(dump);
  fflush(stdout);
}
static void power(arb_t e, int h, int k, int j) {
  arb_mul_si(e, up.p0, k, p);
  arb_addmul_ui(e, up.p1, j, p);
  arb_t t;
  arb_init(t);
  arb_set_si(t, h);
  arb_mul_2exp_si(t, t, -1);
  arb_add(e, e, t, p);
  arb_clear(t);
}
/* A smooth Fourier correction sum with an enclosing infinite tail. */
static void smooth_sum(arb_t out, const arb_t exponent) {
  arb_t n, t, d;
  arb_init(n);
  arb_init(t);
  arb_init(d);
  arb_zero(out);
  for (int m = 1; m <= TANH_MAX; m++) {
    arb_set_ui(n, m);
    arb_tanh(d, n, p);
    arb_sqrt(d, d, p);
    arb_neg(d, d);
    arb_add_ui(d, d, 1, p);
    arb_pow(t, n, exponent, p);
    arb_addmul(out, d, t, p);
  }
  compute_tail_error(t, exponent, TANH_MAX + 1, p);
  arb_add_error(out, t);
  arb_clear(n);
  arb_clear(t);
  arb_clear(d);
}
/* Full x^(2m) coefficient of Q, including smooth correction and finite b. */
static void qcoef(arb_t out, int m) {
  arb_t q, a, t, e, fac;
  arb_init(q);
  arb_init(a);
  arb_init(t);
  arb_init(e);
  arb_init(fac);
  compute_A_m_0(out, m, &up, p);
  arb_fac_ui(fac, 2 * m, p);
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? up.N1 : up.N0); k++) {
      power(q, 0, k, j);
      arb_neg(e, q);
      arb_add_si(e, e, 2 * m - 2, p);
      smooth_sum(t, e);
      arb_mul(t, t, (j ? up.a1k : up.a0k) + k, p);
      arb_div(t, t, fac, p);
      /* d=sqrt(tanh)-1; coefficient (1-cos) */
      if (m % 2)
        arb_neg(t, t);
      arb_add(out, out, t, p);
    }
  for (int n = 1; n <= up.N2; n++) {
    arb_set_ui(t, n);
    arb_tanh(a, t, p);
    arb_sqrt(a, a, p);
    arb_set_si(e, 4 * m - 1);
    arb_mul_2exp_si(e, e, -1);
    arb_pow(t, t, e, p);
    arb_mul(t, t, a, p);
    arb_mul(t, t, up.bi + n - 1, p);
    arb_div(t, t, fac, p);
    if (m % 2)
      arb_neg(t, t);
    arb_add(out, out, t, p);
  }
  arb_clear(q);
  arb_clear(a);
  arb_clear(t);
  arb_clear(e);
  arb_clear(fac);
}
/* Corrected E_Lu: start with inherited tail, repair n^(3/2)->n^(11/2),
 * add the omitted infinite smooth correction remainder. */
static void elu(arb_t out, const arb_t x) {
  arb_t q, t, a, b, e;
  arb_init(q);
  arb_init(t);
  arb_init(a);
  arb_init(b);
  arb_init(e);
  compute_E_x_Lu0(out, x, 7, &up, p);
  for (int n = 1; n <= up.N2; n++) {
    arb_set_ui(a, n);
    arb_tanh(t, a, p);
    arb_sqrt(t, t, p);
    arb_sqrt(a, a, p);
    arb_mul_ui(a, a, n, p);
    arb_mul(a, a, t, p);
    arb_mul_ui(a, a, (ulong)n * n * n * n - 1, p);
    arb_abs(b, up.bi + n - 1);
    arb_mul(a, a, b, p);
    arb_div_ui(a, a, 720, p);
    arb_pow_ui(b, x, 6, p);
    arb_addmul(out, a, b, p);
  }
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? up.N1 : up.N0); k++) {
      power(q, 0, k, j);
      arb_neg(e, q);
      arb_add_ui(e, e, 4, p);
      smooth_sum(t, e);
      arb_abs(t, t);
      arb_abs(a, (j ? up.a1k : up.a0k) + k);
      arb_mul(t, t, a, p);
      arb_div_ui(t, t, 720, p);
      arb_pow_ui(b, x, 6, p);
      arb_addmul(out, t, b, p);
    }
  arb_clear(q);
  arb_clear(t);
  arb_clear(a);
  arb_clear(b);
  arb_clear(e);
}
typedef struct {
  int h, k, j;
  arb_t c;
} term;
static term ds[120];
static int nd;
static void addterm(int h, int k, int j, const arb_t a) {
  int i;
  for (i = 0; i < nd; i++)
    if (ds[i].h == h && ds[i].k == k && ds[i].j == j)
      break;
  if (i == nd) {
    nd++;
    ds[i].h = h;
    ds[i].k = k;
    ds[i].j = j;
    arb_init(ds[i].c);
  }
  arb_add(ds[i].c, ds[i].c, a, p);
}
static int cusp(arb_t out, const arb_t eps) {
  term ps[10];
  arb_t a, e, t, pu, eu, el, lo, lambda;
  arb_init(a);
  arb_init(e);
  arb_init(t);
  arb_init(pu);
  arb_init(eu);
  arb_init(el);
  arb_init(lo);
  arb_init(lambda);
  int z = 0;
  for (int j = 0; j < 2; j++)
    for (int k = 0; k < (j ? up.N1 : up.N0); k++) {
      ps[z].h = 1;
      ps[z].k = k;
      ps[z].j = j;
      arb_init(ps[z].c);
      compute_a_jk_0(ps[z].c, j, k, &up, p);
      compute_A_jk_0(a, j, k, &up, p);
      addterm(2, k, j, a);
      z++;
    }
  for (int m = 1; m <= 2; m++) {
    ps[z].h = 4 * m;
    ps[z].k = ps[z].j = 0;
    arb_init(ps[z].c);
    if (m == 1)
      compute_b_0_0(ps[z].c, &up, p);
    else
      compute_b_1_0(ps[z].c, &up, p);
    qcoef(a, m);
    addterm(4 * m, 0, 0, a);
    z++;
  }
  for (int i = 0; i < z; i++)
    for (int j = 0; j < z; j++) {
      arb_mul(a, ps[i].c, ps[j].c, p);
      arb_neg(a, a);
      addterm(ps[i].h + ps[j].h, ps[i].k + ps[j].k, ps[i].j + ps[j].j, a);
    }
  arb_zero(out);
  arb_zero(pu);
  arb_const_pi(lambda, p);
  arb_div_ui(lambda, lambda, 8, p);
  arb_sqrt(lambda, lambda, p);
  arb_set(lo, lambda);
  int ok = 1;
  for (int i = 0; i < nd; i++) {
    if (ds[i].h == 2 && ds[i].k == 0 && ds[i].j == 0) {
      /* Exact identity pi*a00/2=2*pi*a00^2 when a00=1/4. */
      arb_set_str(a, "0.25", p);
      if (!arb_equal(a, up.a0k))
        ok = 0;
      continue;
    }
    power(e, ds[i].h - 3, ds[i].k, ds[i].j);
    if (!arb_is_positive(e))
      ok = 0;
    arb_pow(t, eps, e, p);
    arb_abs(a, ds[i].c);
    arb_addmul(out, a, t, p);
  }
  for (int i = 0; i < z; i++) {
    power(e, ps[i].h - 1, ps[i].k, ps[i].j);
    arb_pow(t, eps, e, p);
    arb_abs(a, ps[i].c);
    arb_mul(a, a, t, p);
    arb_add(pu, pu, a, p);
    if (i > 0)
      arb_sub(lo, lo, a, p);
  }
  compute_E_x_u0(eu, eps, 7, &up, p);
  elu(el, eps);
  arb_sqrt(t, eps, p);
  arb_div(a, eu, t, p);
  arb_sub(lo, lo, a, p);
  /* Errors divided by x^(3/2); each majorant is a sum of positive powers,
    hence its supremum after scaling occurs at eps. */
  arb_mul(t, t, eps, p);
  arb_div(el, el, t, p);
  arb_add(out, out, el, p);
  arb_mul(a, pu, eu, p);
  arb_div(a, a, eps, p);
  arb_mul_ui(a, a, 2, p);
  arb_add(out, out, a, p);
  arb_sqr(a, eu, p);
  arb_div(a, a, t, p);
  arb_add(out, out, a, p);
  arb_div(out, out, lo, p);
  arb_mul_2exp_si(out, out, -1);
  show("cusp lower u0/sqrt(x)", lo);
  show("cusp |F0| bound", out);
  ok = ok && arb_is_positive(lo) && arb_is_finite(out);
  /* Appendix B prints lambda_- > 0.604 and the cusp defect < 1.55e-8. */
  arb_set_str(a, "0.604", p);
  ok = ok && arb_gt(lo, a);
  arb_set_str(a, "1.55e-8", p);
  ok = ok && arb_lt(out, a);
  for (int i = 0; i < z; i++)
    arb_clear(ps[i].c);
  for (int i = 0; i < nd; i++)
    arb_clear(ds[i].c);
  nd = 0;
  arb_clear(a);
  arb_clear(e);
  arb_clear(t);
  arb_clear(pu);
  arb_clear(eu);
  arb_clear(el);
  arb_clear(lo);
  arb_clear(lambda);
  return ok;
}
static long leaves = 0, nodes = 0;
static int compact(const arb_t l, const arb_t r, const arb_t target, int depth) {
  arb_t x, m, f, d, t, al, u, rad;
  arb_init(x);
  arb_init(m);
  arb_init(f);
  arb_init(d);
  arb_init(t);
  arb_init(al);
  arb_init(u);
  arb_init(rad);
  arb_union(x, l, r, p);
  arb_add(m, l, r, p);
  arb_mul_2exp_si(m, m, -1);
  arb_sub(rad, x, m, p);
  compute_F0(f, m, &up, p);
  compute_dF0(d, x, &up, p);
  arb_mul(t, d, rad, p);
  arb_add(f, f, t, p);
  arb_abs(f, f);
  compute_u0(u, x, &up, p);
  arb_mul_ui(t, u, 2, p);
  arb_div(al, x, t, p);
  arb_set_str(t, "2.7", p);
  int ok = arb_is_positive(u) && arb_lt(f, target) && arb_lt(al, t);
  nodes++;
  if (ok) {
    leaves++;
    if (leaves % 1000 == 0) {
      printf("compact leaves=%ld nodes=%ld x=", leaves, nodes);
      arb_printn(m, 8, 0);
      puts("");
      fflush(stdout);
    }
  } else if (depth < 40 && nodes < 2000000) {
    ok = compact(l, m, target, depth + 1);
    if (ok)
      ok = compact(m, r, target, depth + 1);
  } else {
    show("UNRESOLVED cell", x);
    show("F0 enclosure", f);
    show("alpha enclosure", al);
  }
  arb_clear(x);
  arb_clear(m);
  arb_clear(f);
  arb_clear(d);
  arb_clear(t);
  arb_clear(al);
  arb_clear(u);
  arb_clear(rad);
  return ok;
}
int main(int argc, char **argv) {
  populate_parameters(&up, p);
  arb_t eps, d, pi, v, t;
  arb_init(eps);
  arb_init(d);
  arb_init(pi);
  arb_init(v);
  arb_init(t);
  arb_const_pi(pi, p);
  arb_set_str(eps, argc > 1 ? argv[1] : "0.001", p);
  arb_set_str(d, "0.0000001", p);
  puts("Existence u0/F0 verifier: inherited exact-decimal input variant");
  if (argc > 3 && !strcmp(argv[3], "--reject-test"))
    arb_zero(d);
  int ok = cusp(v, eps) && arb_lt(v, d);
  compute_validation_T0(v, eps, 7, &up, p);
  show("local T0 scalar (CB replaced by 1)", v);
  if (argc > 2) {
    ok = ok && compact(eps, pi, d, 0);
    printf("u0/F0 compact coverage: %s leaves=%ld nodes=%ld\n", ok ? "PASS" : "UNRESOLVED", leaves,
           nodes);
  } else
    printf("cusp defect target: %s\n", ok ? "PASS" : "UNRESOLVED");
  puts("This executable certifies only u0/F0; run the complete certificate entrypoint for "
       "existence.");
  return ok ? 0 : 2;
}
