/* norm_t0.c: the operator norm of T_0 on the compact cells (stage T0).
 *
 * Evaluates the row norm of T_0 through its primitive on a recursive cover
 * of [eps, pi] and certifies the bound 0.998 of the proposition on T_0 cell
 * by cell, printing the cover record checked by run.py. The cusp interval is
 * handled by the scalar stages through the local criterion (ineqT0).
 */
#include "constant_functions.h"
#include "complex_functions.h"
#include <stdlib.h>
static slong p = CAP_PREC;
static u_params_struct up;
static long nodes, leaves;
static arb_t norm_limit;
static void show(const char *s, const arb_t x) {
  printf("%s = ", s);
  arb_printn(x, 22, 0);
  puts("");
  fflush(stdout);
}
static void g(arb_t out, const arb_t x, const arb_t w) {
  arb_t y, t, a;
  arb_init(y);
  arb_init(t);
  arb_init(a);
  arb_mul(y, x, w, p);
  arb_sub(t, x, y, p);
  Kp(out, t, p);
  arb_add(t, x, y, p);
  Kp(a, t, p);
  arb_add(out, out, a, p);
  Kp(a, y, p);
  arb_mul_ui(a, a, 2, p);
  arb_sub(out, out, a, p);
  arb_clear(y);
  arb_clear(t);
  arb_clear(a);
}
/* J_x(y) with dJ/dy = y G_x(y). */
static void primitive(arb_t out, const arb_t x, const arb_t y) {
  arb_t t, a, b;
  arb_init(t);
  arb_init(a);
  arb_init(b);
  arb_sub(t, x, y, p);
  IKp(a, t, p);
  arb_mul(out, y, a, p);
  IIKp(a, t, p);
  arb_add(out, out, a, p);
  arb_neg(out, out);
  arb_add(t, x, y, p);
  IKp(a, t, p);
  arb_addmul(out, y, a, p);
  IIKp(a, t, p);
  arb_sub(out, out, a, p);
  IKp(a, y, p);
  arb_mul_ui(a, a, 2, p);
  arb_submul(out, y, a, p);
  IIKp(a, y, p);
  arb_mul_ui(a, a, 2, p);
  arb_add(out, out, a, p);
  arb_clear(t);
  arb_clear(a);
  arb_clear(b);
}
static int row(arb_t out, const arb_t x) {
  arb_t l, r, m, v, z, t, a, pi, u;
  arb_init(l);
  arb_init(r);
  arb_init(m);
  arb_init(v);
  arb_init(z);
  arb_init(t);
  arb_init(a);
  arb_init(pi);
  arb_init(u);
  arb_const_pi(pi, p);
  arb_set_str(l, "0.4", p);
  arb_set_str(r, "0.95", p);
  g(v, x, l);
  int ok = arb_is_negative(v);
  g(v, x, r);
  ok = ok && arb_is_positive(v);
  for (int k = 0; k < 14 && ok; k++) {
    arb_add(m, l, r, p);
    arb_mul_2exp_si(m, m, -1);
    g(v, x, m);
    if (arb_is_negative(v))
      arb_set(l, m);
    else if (arb_is_positive(v))
      arb_set(r, m);
    else {
      /* Mean value theorem and G_w >= -x K'(x*r) > 0 on [l,r].
    This encloses all roots as x ranges over the cell, even when a
    bisection midpoint meets the root. */
      arb_mul(t, x, r, p);
      dKp(a, t, p);
      arb_neg(a, a);
      arb_mul(a, a, x, p);
      if (arb_is_positive(a)) {
        arb_abs(v, v);
        arb_div(v, v, a, p);
        arb_set(z, m);
        arb_add_error(z, v);
        arb_union(t, l, r, p);
        if (!arb_intersection(z, z, t, p)) {
          ok = 0;
          break;
        }
        arb_get_lbound_arf(arb_midref(l), z, p);
        mag_zero(arb_radref(l));
        arb_get_ubound_arf(arb_midref(r), z, p);
        mag_zero(arb_radref(r));
      }
      break;
    }
  }
  arb_union(z, l, r, p);
  arb_mul(z, z, x, p);
  primitive(out, x, z);
  arb_mul_si(out, out, -2, p);
  IIKp(a, x, p);
  arb_mul_ui(a, a, 2, p);
  arb_sub(out, out, a, p);
  arb_sub(t, pi, x, p);
  arb_abs(t, t);
  IIKp(a, t, p);
  arb_mul_ui(a, a, 2, p);
  arb_sub(out, out, a, p);
  arb_sub(out, out, x, p);
  IIKp(a, pi, p);
  arb_mul_ui(a, a, 2, p);
  arb_add(out, out, a, p);
  compute_u0(u, x, &up, p);
  arb_mul(t, x, u, p);
  arb_mul_ui(t, t, 2, p);
  arb_div(out, out, t, p);
  ok = ok && arb_is_positive(u) && arb_is_finite(out);
  arb_clear(l);
  arb_clear(r);
  arb_clear(m);
  arb_clear(v);
  arb_clear(z);
  arb_clear(t);
  arb_clear(a);
  arb_clear(pi);
  arb_clear(u);
  return ok;
}
static int cover(const arb_t l, const arb_t r, int depth) {
  arb_t x, v, m;
  arb_init(x);
  arb_init(v);
  arb_init(m);
  arb_union(x, l, r, p);
  nodes++;
  int ok = row(v, x) && arb_lt(v, norm_limit);
  if (ok) {
    leaves++;
    if (leaves % 100 == 0) {
      printf("T0 leaves=%ld nodes=%ld x=", leaves, nodes);
      arb_printn(x, 8, 0);
      puts("");
      fflush(stdout);
    }
  } else if (depth < 40 && nodes < 1000000) {
    arb_add(m, l, r, p);
    arb_mul_2exp_si(m, m, -1);
    ok = cover(l, m, depth + 1);
    if (ok)
      ok = cover(m, r, depth + 1);
  } else {
    show("UNRESOLVED x", x);
    show("row", v);
  }
  arb_clear(x);
  arb_clear(v);
  arb_clear(m);
  return ok;
}
int main(int argc, char **argv) {
  populate_parameters(&up, p);
  arb_init(norm_limit);
  arb_set_str(norm_limit, "0.998", p);
  arb_t eps, pi, t;
  arb_init(eps);
  arb_init(pi);
  arb_init(t);
  arb_set_str(eps, "0.05", p);
  arb_const_pi(pi, p);
  if (argc > 1) {
    arb_set_str(eps, argv[1], p);
    int ok = row(t, eps);
    show("ROW DIAGNOSTIC", t);
    return ok ? 0 : 2;
  }
  int ok = cover(eps, pi, 0);
  printf("T0 compact coverage %s leaves=%ld nodes=%ld\n", ok ? "PASS" : "UNRESOLVED", leaves,
         nodes);
  puts("Local T0 and scalar assembly are separate obligations");
  return ok ? 0 : 2;
}
