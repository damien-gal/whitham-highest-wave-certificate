/* test_kernel.c: regression output of the kernel evaluators (stage kernel_tests).
 * These point values are compared with an independent implementation; they are
 * not part of the proof.
 */
#include "kernel.h"
int main(void) {
  Kernel k;
  kernel_init(&k, 80, 40);
  arb_ptr cs[2] = {_arb_vec_init(13), _arb_vec_init(13)};
  root_coefficients(cs[0], cs[1]);
  const char *points[] = {".01", ".1", "1", "2", "3.14"};
  for (int i = 0; i < 5; i++) {
    size_t m = mark();
    B x = decimal(points[i]);
    for (int s = 0; s < 2; s++) {
      arb_ptr v = tmp(), d = tmp(), c = tmp();
      poly_box(v, d, c, cs[s], 13, x, x);
      printf("POINT %s %d\n", points[i], s);
      show("W", curve(&k, x, v));
      show("Wprime", curve_derivative(&k, x, v, d));
    }
    reset(m);
  }
  show("M", regular_second(&k, add(k.pi, decimal(".001")), 0));
  show("M3", regular_second(&k, times(k.pi, 2), 1));
  require(lattice(rational(-1, 3), 4, 0) == -6 && lattice(rational(-1, 3), 4, 1) == -5,
          "negative lattice rounding");
  require(lattice(rational(1, 3), 4, 0) == 5 && lattice(rational(1, 3), 4, 1) == 6,
          "positive lattice rounding");
  puts("PASS ball conversion tests");
  return 0;
}
