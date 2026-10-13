/* matrix_driver.c: kernel moment matrices on the refinement grid (stage matrix).
 *
 * Usage: matrix <grid.bin> <N> <output.bin>. The grid is N increasing
 * binary64 numbers whose last entry is the double nearest to pi; the output
 * holds outward binary64 enclosures of the kernel moments and of the seed
 * values of u_0 on the grid. The moments are computed with the kernel
 * primitives of core/basic_functions.c and rounded outward; they are the
 * input of uniqueness.c and are not a certificate by themselves.
 */
#include "basic_functions.h"
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static slong prec = CAP_PREC;
static void die(const char *s) {
  fprintf(stderr, "FAIL: %s\n", s);
  exit(2);
}
static void kernel(arb_t v, const arb_t x, const arb_t z, const arb_t y) {
  arb_t t, a;
  arb_init(t);
  arb_init(a);
  arb_sub(t, x, y, prec);
  Kp(v, t, prec);
  arb_add(t, x, y, prec);
  Kp(a, t, prec);
  arb_add(v, v, a, prec);
  arb_sub(t, z, y, prec);
  Kp(a, t, prec);
  arb_sub(v, v, a, prec);
  arb_add(t, z, y, prec);
  Kp(a, t, prec);
  arb_sub(v, v, a, prec);
  arb_clear(t);
  arb_clear(a);
}
static void integral(arb_t out, const arb_t x) {
  if (arb_is_zero(x)) {
    arb_zero(out);
    return;
  }
  if (arb_contains_zero(x)) {
    arb_t ax;
    arb_init(ax);
    arb_abs(ax, x);
    arb_get_ubound_arf(arb_midref(ax), ax, prec);
    mag_zero(arb_radref(ax));
    arb_t val;
    arb_init(val);
    IKp(val, ax, prec);
    arb_zero(out);
    arb_add_error(out, val);
    arb_clear(val);
    arb_clear(ax);
    return;
  }
  if (arb_is_negative(x)) {
    arb_t a;
    arb_init(a);
    arb_neg(a, x);
    integral(out, a);
    arb_neg(out, out);
    arb_clear(a);
    return;
  }
  arb_t pi, a;
  arb_init(pi);
  arb_init(a);
  arb_const_pi(pi, prec);
  if (arb_gt(x, pi)) {
    arb_mul_2exp_si(a, pi, 1);
    arb_sub(a, x, a, prec);
    integral(out, a);
    arb_add_ui(out, out, 1, prec);
    arb_clear(pi);
    arb_clear(a);
    return;
  }
  arb_clear(pi);
  arb_clear(a);
  IKp(out, x, prec);
  if (!arb_is_finite(out)) {
    fprintf(stderr, "IKp argument: ");
    arb_printn(x, 20, 0);
    die("IKp nonfinite");
  }
}
static void J(arb_t v, const arb_t x, const arb_t y) {
  arb_t t, a;
  arb_init(t);
  arb_init(a);
  arb_sub(t, x, y, prec);
  integral(v, t);
  arb_neg(v, v);
  arb_add(t, x, y, prec);
  integral(a, t);
  arb_add(v, v, a, prec);
  integral(a, y);
  arb_mul_ui(a, a, 2, prec);
  arb_sub(v, v, a, prec);
  arb_clear(t);
  arb_clear(a);
}
static void prim(arb_t v, const arb_t x, const arb_t z, const arb_t y) {
  arb_t a;
  arb_init(a);
  J(v, x, y);
  if (!arb_is_zero(z)) {
    J(a, z, y);
    arb_sub(v, v, a, prec);
  }
  arb_clear(a);
}
static void emit(FILE *f, const arb_t v, int nonnegative) {
  if (!arb_is_finite(v))
    die("nonfinite output");
  arf_t lo, hi;
  arf_init(lo);
  arf_init(hi);
  arb_get_interval_arf(lo, hi, v, prec);
  double d[2] = {arf_get_d(lo, ARF_RND_FLOOR), arf_get_d(hi, ARF_RND_CEIL)};
  if (nonnegative) {
    d[0] = fmax(d[0], 0);
    d[1] = fmax(d[1], 0);
  }
  if (!isfinite(d[0]) || !isfinite(d[1]) || d[0] > d[1])
    die("invalid endpoints");
  if (fwrite(d, sizeof(double), 2, f) != 2)
    die("write");
  arf_clear(lo);
  arf_clear(hi);
}
static void row(FILE *f, arb_srcptr grid, int n, int i, int k) {
  const arb_struct *x = grid + i + 1, *z = grid + k + 1;
  arb_t lo, hi, mid, t, v, root, jroot, left, right, pos, neg;
  arb_init(lo);
  arb_init(hi);
  arb_init(mid);
  arb_init(t);
  arb_init(v);
  arb_init(root);
  arb_init(jroot);
  arb_init(left);
  arb_init(right);
  arb_init(pos);
  arb_init(neg);
  arb_sub(t, x, z, prec);
  arb_mul_2exp_si(t, t, -16);
  arb_add(lo, z, t, prec);
  arb_sub(hi, x, t, prec);
  kernel(v, x, z, lo);
  if (!arb_is_negative(v))
    die("root left sign");
  kernel(v, x, z, hi);
  if (!arb_is_positive(v))
    die("root right sign");
  for (int b = 0; b < 55; b++) {
    arb_add(mid, lo, hi, prec);
    arb_mul_2exp_si(mid, mid, -1);
    kernel(v, x, z, mid);
    if (arb_is_negative(v))
      arb_set(lo, mid);
    else if (arb_is_positive(v))
      arb_set(hi, mid);
    else {
      arb_sub(t, x, lo, prec);
      dKp(t, t, prec);
      arb_neg(t, t);
      if (!arb_is_positive(t))
        die("root derivative bound");
      arb_abs(v, v);
      arb_div(v, v, t, prec);
      arb_set(root, mid);
      arb_add_error(root, v);
      arb_union(t, lo, hi, prec);
      if (!arb_intersection(root, root, t, prec))
        die("root intersection");
      arb_get_lbound_arf(arb_midref(lo), root, prec);
      mag_zero(arb_radref(lo));
      arb_get_ubound_arf(arb_midref(hi), root, prec);
      mag_zero(arb_radref(hi));
      break;
    }
  }
  arb_union(root, lo, hi, prec);
  prim(jroot, x, z, root);
  arb_zero(left);
  for (int j = 0; j < n; j++) {
    if (j == n - 1)
      arb_zero(right);
    else
      prim(right, x, z, grid + j + 1);
    if (arb_le(grid + j + 1, lo)) {
      arb_zero(pos);
      arb_sub(neg, left, right, prec);
    } else if (arb_ge(grid + j, hi)) {
      arb_sub(pos, right, left, prec);
      arb_zero(neg);
    } else {
      if (!arb_le(grid + j, lo) || !arb_ge(grid + j + 1, hi))
        die("root overlaps a grid boundary; subdivide bracket");
      arb_sub(pos, right, jroot, prec);
      arb_sub(neg, left, jroot, prec);
    }
    emit(f, pos, 1);
    emit(f, neg, 1);
    arb_set(left, right);
  }
  arb_clear(lo);
  arb_clear(hi);
  arb_clear(mid);
  arb_clear(t);
  arb_clear(v);
  arb_clear(root);
  arb_clear(jroot);
  arb_clear(left);
  arb_clear(right);
  arb_clear(pos);
  arb_clear(neg);
}
int main(int argc, char **argv) {
  if (argc != 4 && argc != 5)
    die("usage grid N output");
  int n = atoi(argv[2]);
  if (n < 2 || n > 10000)
    die("N");
  FILE *f = fopen(argv[1], "rb");
  if (!f)
    die("grid open");
  double *d = malloc(n * sizeof(double));
  if (fread(d, sizeof(double), n, f) != (size_t)n)
    die("grid read");
  fclose(f);
  arb_ptr grid = _arb_vec_init(n + 1);
  for (int i = 0; i < n; i++) {
    if (!(d[i] > 0) || !isfinite(d[i]) || (i && d[i] <= d[i - 1]))
      die("grid order");
    arb_set_d(grid + i + 1, d[i]);
  }
  arb_const_pi(grid + n, prec);
  f = fopen(argv[3], "wb");
  if (!f)
    die("output open");
  for (int i = 0; i < n; i++) {
    if (argc == 5 && (strcmp(argv[4], "seeds") == 0 || i != atoi(argv[4])))
      continue;
    row(f, grid, n, i, -1);
    int k = 0;
    while (k + 1 < i && d[k + 1] < .4 * d[i])
      k++;
    if (i)
      row(f, grid, n, i, k);
    else {
      double z = 0;
      for (int j = 0; j < 4 * n; j++)
        fwrite(&z, 8, 1, f);
    }
    if (i % 20 == 0) {
      fprintf(stderr, "row %d / %d\n", i, n);
      fflush(f);
    }
  }
  if (argc == 5 && strcmp(argv[4], "seeds") != 0) {
    fclose(f);
    return 0;
  }
  arb_t z, a, b, h;
  arb_init(z);
  arb_init(a);
  arb_init(b);
  arb_init(h);
  u_params_struct up;
  populate_parameters(&up, prec);
  for (int i = 0; i < n; i++) {
    arb_set_str(z, "0.437", prec);
    arb_mul(z, z, grid + i + 1, prec);
    prim(a, grid + i + 1, z, z);
    arb_neg(a, a);
    prim(b, grid + i + 1, z, grid + i + 1);
    arb_neg(b, b);
    arb_min(h, a, b, prec);
    arb_mul_2exp_si(h, h, -1);
    emit(f, h, 1);
    compute_u0(a, grid + i + 1, &up, prec);
    emit(f, a, 0);
  }
  if (fclose(f))
    die("close");
  fprintf(stderr, "PASS: kernel moment and seed enclosures only\n");
  return 0;
}
