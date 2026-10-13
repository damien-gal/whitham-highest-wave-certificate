/* test_evaluators.c: independent cross-checks of the special-function
 * evaluators (stage evaluators). Numerical agreement with an independent
 * evaluation is a regression test, not the proof of the evaluators; the
 * rejection paths are exercised as well.
 */
#include "basic_functions.h"
int main(void) {
  const slong p = CAP_PREC;
  const char *orders[] = {"-1.5", "-0.5", "0.5", "1.5", "2.5", "3.2", "4.87360289266183309"};
  const char *points[] = {"0.00001",          "0.05", "0.3", "2", "3.141592653589793", "5.8",
                          "6.283185307179586"};
  arb_t s, x, a, b, t, rad, wide;
  arb_init(s);
  arb_init(x);
  arb_init(a);
  arb_init(b);
  arb_init(t);
  arb_init(rad);
  arb_init(wide);
  acb_t z, w;
  acb_init(z);
  acb_init(w);
  int checks = 0, ok = 1;
  for (int i = 0; i < 7; i++)
    for (int j = 0; j < 7; j++)
      for (int sine = 0; sine < 2; sine++) {
        arb_set_str(s, orders[i], p);
        arb_set_str(x, points[j], p);
        ok = ok && clausen_cached(a, s, x, sine, p);
        arb_cos(acb_realref(z), x, p + 32);
        arb_sin(acb_imagref(z), x, p + 32);
        acb_set_arb(w, s);
        acb_polylog(w, w, z, p + 32);
        arb_set(b, sine ? acb_imagref(w) : acb_realref(w));
        if (!arb_is_finite(a) || !arb_is_finite(b) || !arb_overlaps(a, b)) {
          printf("MISMATCH s=%s x=%s sine=%d\n", orders[i], points[j], sine);
          ok = 0;
        }
        /* A full interval enclosure must also contain enclosures at its midpoint. */
        arb_set(wide, x);
        arb_set_str(rad, "0.00000000000000000001", p);
        arb_add_error(wide, rad);
        if (!clausen_cached(t, s, wide, sine, p) || !arb_overlaps(t, b))
          ok = 0;
        checks++;
      }
  /* Exact pi as a ball must not trigger a false domain rejection. */
  arb_const_pi(x, p);
  arb_set_str(s, "2.5", p);
  if (!clausen_cached(a, s, x, 0, p))
    ok = 0;
  arb_set_str(s, "0.5", p);
  arb_set_ui(x, 20);
  if (clausen_cached(a, s, x, 0, p))
    ok = 0;
  arb_set_si(s, -3);
  arb_set_ui(x, 1);
  if (clausen_cached(a, s, x, 0, p))
    ok = 0;
  printf("Clausen direct-polylog/interval regression checks: %s (%d cases)\n", ok ? "PASS" : "FAIL",
         checks);
  return ok ? 0 : 2;
}
