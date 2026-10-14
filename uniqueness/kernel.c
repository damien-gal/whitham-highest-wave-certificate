/* kernel.c: enclosures of the kernel K and of the quantities of Appendix C.
 *
 * regular():        the regular part of K as a series of Clausen coefficients
 *                   zeta(1/2-2j)/(2j)!, with the analytic tail of Appendix C;
 * homogeneous():    the homogeneous cusp part of the normalized kernel;
 * correction():     the Fourier correction with weights (1-sqrt(tanh n))/sqrt(n),
 *                   truncated after N terms with an explicit geometric tail;
 * curve(), curve_derivative(), small_curve(): the normalized kernel along a
 *                   curve y = x t, its derivative and its small-x form;
 * fourier_second(), regular_second(): bounds on second derivatives, the latter
 *                   with the three-image coefficients when images is set;
 * root_coefficients(): the 26 exact decimal coefficients of the polynomials
 *                   y*_-(x), y*_+(x) printed in Appendix C (they are compared
 *                   with the paper and with inputs/kernel-sign.json by run.py);
 * poly_box():       value, derivative and center of a polynomial on a cell,
 *                   through a centered Taylor shift with interval coefficients.
 */
#include "kernel.h"
void kernel_init(Kernel *k, int M, int N) {
  k->M = M;
  k->N = N;
  arb_init(k->pi);
  arb_init(k->C);
  arb_init(k->lambda);
  arb_init(k->A);
  arb_const_pi(k->pi, CAP_PREC);
  size_t m = mark();
  arb_set(k->C, divide(integer(1), root(times(k->pi, 2))));
  arb_set(k->lambda, root(over(k->pi, 8)));
  arb_set(k->A, divide(zetab(rational(5, 2)), root(k->pi)));
  k->coeff = _arb_vec_init(M + 1);
  k->w = _arb_vec_init(N);
  for (int j = 1; j <= M; j++) {
    arb_ptr fac = tmp();
    arb_fac_ui(fac, 2 * j, CAP_PREC);
    arb_set(k->coeff + j, divide(zetab(rational(1 - 4 * j, 2)), fac));
    if (j % 2)
      arb_neg(k->coeff + j, k->coeff + j);
    require(positive(k->coeff + j), "positive Clausen coefficient");
  }
  for (int n = 1; n <= N; n++) {
    arb_ptr t = tmp();
    arb_tanh(t, integer(n), CAP_PREC);
    arb_set(k->w + n - 1, divide(sub(integer(1), root(t)), root(integer(n))));
  }
  reset(m);
}
void kernel_clear(Kernel *k) {
  _arb_vec_clear(k->coeff, k->M + 1);
  _arb_vec_clear(k->w, k->N);
  arb_clear(k->pi);
  arb_clear(k->C);
  arb_clear(k->lambda);
  arb_clear(k->A);
}
arb_ptr regular(Kernel *k, B z, int derivative) {
  B q = divide(absupper(z), times(k->pi, 2)), r = sq(q), d = sub(integer(1), r);
  require(lt(q, integer(1)), "kernel series radius");
  B v = integer(0);
  for (int j = k->M; j >= 1; j--)
    v = add(mul(v, sq(z)), derivative ? times(k->coeff + j, 2 * j) : k->coeff + j);
  v = mul(v, derivative ? z : sq(z));
  B tail = derivative
               ? divide(mul(mul(divide(k->A, times(k->pi, 2)), times(power(q, 2 * k->M + 1), 2)),
                            sub(integer(k->M + 1), times(r, k->M))),
                        sq(d))
               : divide(mul(k->A, power(q, 2 * k->M + 2)), d);
  return inflate(v, tail);
}
arb_ptr homogeneous(Kernel *k, B t) {
  return mul(k->C,
             sub(add(root(divide(t, sub(integer(1), t))), root(divide(t, add(integer(1), t)))),
                 integer(2)));
}
static arb_ptr correction(Kernel *k, B x, B y, B yp, int derivative) {
  B v = integer(0);
  for (int n = 1; n <= k->N; n++) {
    B nx = times(x, n), ny = times(y, n);
    B term = derivative ? times(add(mul(sinb(nx), cosb(ny)),
                                    mul(mul(sub(integer(1), cosb(nx)), sinb(ny)), neg(yp))),
                                n)
                        : mul(sub(integer(1), cosb(nx)), cosb(ny));
    v = add(v, mul(k->w + n - 1, term));
  }
  B rho = expb(integer(-2)), d = sub(integer(1), rho);
  int n = k->N + 1;
  B tail = derivative ? divide(mul(times(add(integer(1), times(ab(yp), 2)), 4),
                                   mul(power(rho, n), sub(integer(n), times(rho, n - 1)))),
                               mul(k->pi, sq(d)))
                      : divide(times(power(rho, n), 8), mul(k->pi, d));
  return inflate(divide(times(v, 2), k->pi), tail);
}
arb_ptr curve(Kernel *k, B x, B t) {
  B y = mul(x, t);
  B R = divide(
      sub(add(regular(k, sub(x, y), 0), regular(k, add(x, y), 0)), times(regular(k, y, 0), 2)),
      k->pi);
  return add(homogeneous(k, t), mul(root(y), add(R, correction(k, x, y, NULL, 0))));
}
arb_ptr curve_derivative(Kernel *k, B x, B t, B tp) {
  B y = mul(x, t), yp = add(t, mul(x, tp));
  B R = divide(
      sub(add(regular(k, sub(x, y), 0), regular(k, add(x, y), 0)), times(regular(k, y, 0), 2)),
      k->pi);
  B Rp = divide(sub(add(mul(sub(integer(1), yp), regular(k, sub(x, y), 1)),
                        mul(add(integer(1), yp), regular(k, add(x, y), 1))),
                    times(mul(yp, regular(k, y, 1)), 2)),
                k->pi);
  B hp = times(
      mul(k->C,
          add(divide(integer(1), mul(root(t), mul(sub(integer(1), t), root(sub(integer(1), t))))),
              divide(integer(1), mul(root(t), mul(add(integer(1), t), root(add(integer(1), t))))))),
      1);
  hp = over(hp, 2);
  return add(
      add(mul(hp, tp), mul(divide(yp, times(root(y), 2)), add(R, correction(k, x, y, NULL, 0)))),
      mul(root(y), add(Rp, correction(k, x, y, yp, 1))));
}
static arb_ptr sum2(B r, int n) {
  B d = sub(integer(1), r);
  return mul(power(r, n), add(add(divide(integer(n * n), d), divide(times(r, 2 * n), sq(d))),
                              divide(mul(r, add(integer(1), r)), power(d, 3))));
}
arb_ptr fourier_second(Kernel *k) {
  B v = integer(0);
  for (int n = 1; n <= k->N; n++)
    v = add(v, times(k->w + n - 1, n * n));
  return add(v, times(sum2(expb(integer(-2)), k->N + 1), 2));
}
arb_ptr small_curve(Kernel *k, B t, B eps) {
  B y = mul(eps, upper(t)), z = mul(eps, add(integer(1), upper(t)));
  B er = divide(times(absupper(regular(k, z, 0)), 4), k->pi),
    fr = divide(mul(sq(eps), fourier_second(k)), k->pi);
  return inflate(homogeneous(k, t), mul(root(y), add(er, fr)));
}
arb_ptr regular_second(Kernel *k, B xmax, int images) {
  B E = integer(0);
  int M = images ? 48 : k->M;
  B L = times(k->pi, 2);
  for (int j = 1; j <= M; j++) {
    B c = divide(k->coeff + j, k->pi);
    if (images) {
      arb_ptr bin = tmp();
      fmpz_t z;
      fmpz_init(z);
      fmpz_bin_uiui(z, 4 * j, 2 * j);
      arb_set_fmpz(bin, z);
      fmpz_clear(z);
      B term = divide(times(mul(k->C, divide(bin, power(integer(4), 2 * j))), 2),
                      mul(power(L, 2 * j), root(L)));
      c = sub(c, term);
      require(positive(c), "three-image coefficient positive");
    }
    E = add(E, mul(times(c, 2 * j * (2 * j - 1)), power(xmax, 2 * j - 2)));
  }
  B scale = images ? times(k->pi, 4) : L, r = sq(divide(xmax, scale)), d = sub(integer(1), r);
  int n = M + 1;
  B s1 = mul(power(r, n), add(divide(integer(n), d), divide(r, sq(d))));
  B s2 = sum2(r, n);
  B tail = divide(mul(images ? divide(integer(3), k->pi) : divide(k->A, k->pi),
                      sub(times(s2, 4), times(s1, 2))),
                  mul(sq(scale), r));
  return add(add(E, tail), divide(fourier_second(k), k->pi));
}
void poly_box(arb_t value, arb_t deriv, arb_t center, arb_srcptr cs, int n, B a, B b) {
  size_t m = mark();
  B mid = over(add(a, b), 2);
  B h = inflate(integer(0), upper(over(sub(b, a), 2)));
  arb_ptr shifted = _arb_vec_init(n);
  for (int j = 0; j < n; j++) {
    B v = integer(0);
    for (int k = j; k < n; k++) {
      ulong bin = 1;
      for (int t = 1; t <= j; t++)
        bin = bin * (k - t + 1) / t;
      v = add(v, times(mul(cs + k, power(mid, k - j)), bin));
    }
    arb_set(shifted + j, v);
  }
  B v = shifted + n - 1;
  for (int j = n - 2; j >= 0; j--)
    v = add(mul(v, h), shifted + j);
  arb_set(value, v);
  v = integer(0);
  for (int j = n - 1; j >= 1; j--)
    v = add(mul(v, h), times(shifted + j, j));
  arb_set(deriv, v);
  arb_set(center, shifted);
  _arb_vec_clear(shifted, n);
  reset(m);
}

void root_coefficients(arb_ptr minus, arb_ptr plus) {
  require(!arb_set_str(minus + 0, "0.652894520981", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 1, "0.000310430727311", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 2, "-0.0117472533516", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 3, "-0.0491604772863", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 4, "0.02590961183200", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 5, "0.0671953842890", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 6, "-0.1227887419505", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 7, "0.103442958737", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 8, "-0.0536206421817", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 9, "0.0179569474620", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 10, "-0.003787859708335", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 11, "0.0004564193137939", CAP_PREC), "root coefficient");
  require(!arb_set_str(minus + 12, "-0.00002384037506397", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 0, "0.652904520986", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 1, "0.000310430363916", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 2, "-0.01174724776403", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 3, "-0.04916051076109", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 4, "0.02590971469900", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 5, "0.067195198662", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 6, "-0.1227885301853", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 7, "0.103442799936", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 8, "-0.0536205627930", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 9, "0.01795692125156", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 10, "-0.00378785421981", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 11, "0.000456418653356", CAP_PREC), "root coefficient");
  require(!arb_set_str(plus + 12, "-0.00002384034028113", CAP_PREC), "root coefficient");
}
