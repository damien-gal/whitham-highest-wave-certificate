/* uniqueness.c: the a priori enclosures and the comparison of Sections 8-9
 * and Appendix C (stages uniqueness and reject_uniqueness).
 *
 * Usage: uniqueness <grid.bin> <moments.bin> [reject]. Reads the exact
 * binary64 grid and the outward moment enclosures written by matrix_driver.c,
 * checks their sizes, order and conservation, replaces the first-cell moments
 * by weighted bounds, verifies the universal seeds (the constant A of
 * Section 8 below 0.855), runs 240 coupled cusp/nodal refinement steps with
 * the exact integer arithmetic of refine.c, and validates the cusp comparison
 * and every one of the 719 compact cells, printing one COMPARISON_CELL line
 * per cell, one BOUND line per node and the comparison constant Q < 0.982.
 * The same grid-index selection is used here and in matrix_driver.c; no
 * binary double comparison decides a mathematical inequality. With "reject"
 * the comparison target is made impossible and the program must return
 * UNRESOLVED.
 */
#include "kernel.h"
#include <math.h>
#include <string.h>
#include <float.h>
enum { N = 720 };
#define BS_BITS 32
#define CS_BITS 42
#define BS ((int64_t)1 << BS_BITS)
typedef int64_t V;
int advance(int, const V *, const V *, const V *, const V *, const V *, const V *, const int *, V *,
            V *, V, V, int);
static Kernel k;
static arb_ptr xs;
static arb_t delta, rh, aG, aF;
static V *fl, *fh, *gl, *gh, *dpl, *dnh, L[N], U[N];
static int zi[N];
static double grid[N], *raw;
static V maxv(V a, V b) {
  return a > b ? a : b;
}
static V minv(V a, V b) {
  return a < b ? a : b;
}
static B bi(V n) {
  return dyadic(n, BS_BITS);
}
static B ci(V n) {
  return dyadic(n, CS_BITS);
}
static B lowS(B a, B y) {
  return sub(mul(a, asinb(root(divide(y, a)))), root(nonneg(mul(y, sub(a, y)))));
}
static B highS(B a, B y) {
  return add(add(root(nonneg(mul(y, sub(y, a)))), mul(a, acoshb(root(divide(y, a))))),
             over(mul(a, k.pi), 2));
}
static B S(B a, B y) {
  if (arb_is_zero(y))
    return integer(0);
  B pm;
  if (le(y, a))
    pm = lowS(a, y);
  else if (le(a, y))
    pm = highS(a, y);
  else
    pm = hull(lowS(a, minimum_b(y, a)), highS(a, maximum_b(y, a)));
  return mul(k.C, sub(add(pm, root(mul(y, add(y, a)))),
                      add(mul(a, asinhb(root(divide(y, a)))), times(y, 2))));
}
static B ip(B x) {
  if (arb_is_nonnegative(x))
    return times(mul(k.C, root(nonneg(x))), 2);
  if (arb_is_nonpositive(x))
    return neg(times(mul(k.C, root(nonneg(neg(x)))), 2));
  return inflate(integer(0), times(mul(k.C, root(absupper(x))), 2));
}
static B hs(B t, B z) {
  return add(sub(sub(ip(sub(integer(1), t)), ip(add(integer(1), t))), ip(sub(z, t))),
             ip(add(z, t)));
}
static B hsign(B t) {
  return sub(add(divide(integer(1), root(sub(integer(1), t))),
                 divide(integer(1), root(add(integer(1), t)))),
             divide(integer(2), root(t)));
}
static void init_homogeneous(void) {
  arb_init(delta);
  arb_init(rh);
  arb_init(aG);
  arb_init(aF);
  size_t m = mark();
  arb_set(delta, rational(2, 5));
  arb_t a, b;
  arb_init(a);
  arb_init(b);
  arb_set(a, lower(decimal(".65")));
  arb_set(b, upper(decimal(".66")));
  require(negative(hsign(a)) && positive(hsign(b)), "homogeneous root endpoints");
  for (int j = 0; j < 100; j++) {
    size_t q = mark();
    B mid = over(add(a, b), 2), v = hsign(mid);
    if (negative(v))
      arb_set(a, mid);
    else if (positive(v))
      arb_set(b, mid);
    else {
      reset(q);
      break;
    }
    reset(q);
  }
  arb_set(rh, hull(a, b));
  arb_set(aG, neg(S(integer(1), rh)));
  arb_set(aF, add(aG, k.lambda));
  arb_clear(a);
  arb_clear(b);
  reset(m);
}
static void monotone(void) {
  for (int i = 1; i < N; i++)
    L[i] = maxv(L[i], L[i - 1]);
  for (int i = N - 2; i >= 0; i--)
    U[i] = minv(U[i], U[i + 1]);
  for (int i = 0; i < N; i++)
    require(L[i] <= U[i], "ordered profile bounds");
}
static void inputs(const char *gridpath, const char *momentpath) {
  uint16_t endian = 1;
  require(sizeof(double) == 8 && DBL_MANT_DIG == 53 && FLT_RADIX == 2 &&
              *(unsigned char *)&endian == 1,
          "little-endian IEEE binary64 input format");
  FILE *f = fopen(gridpath, "rb");
  require(f != NULL, "grid file");
  require(fread(grid, sizeof(double), N, f) == N && fgetc(f) == EOF, "grid size");
  fclose(f);
  xs = _arb_vec_init(N);
  for (int i = 0; i < N; i++) {
    require(isfinite(grid[i]) && grid[i] > 0 && (!i || grid[i] > grid[i - 1]), "grid order");
    arb_set_d(xs + i, grid[i]);
  }
  arb_set(xs + N - 1, k.pi);
  require(lt(xs + N - 2, k.pi), "last interval ends at pi");
  size_t count = (size_t)N * N * 8 + N * 4;
  raw = malloc(count * sizeof(double));
  require(raw != NULL, "matrix allocation");
  f = fopen(momentpath, "rb");
  require(f != NULL, "moment file");
  require(fread(raw, sizeof(double), count, f) == count && fgetc(f) == EOF, "moment file size");
  fclose(f);
  V **mats[] = {&fl, &fh, &gl, &gh, &dpl, &dnh};
  for (int j = 0; j < 6; j++) {
    *mats[j] = calloc((size_t)N * N, sizeof(V));
    require(*mats[j] != NULL, "integer matrix allocation");
  }
  for (int i = 0; i < N; i++) {
    zi[i] = 0;
    while (zi[i] + 1 < i && grid[zi[i] + 1] < .4 * grid[i])
      zi[i]++;
    for (int mode = 0; mode < 2; mode++) {
      __int128 clo = 0, chi = 0;
      for (int j = 0; j < N; j++) {
        size_t m = mark(), t = ((size_t)i * 2 + mode) * N * 4 + j * 4;
        V pl = lattice(exact_double(raw[t]), 42, 0), ph = lattice(exact_double(raw[t + 1]), 42, 1),
          nl = lattice(exact_double(raw[t + 2]), 42, 0),
          nh = lattice(exact_double(raw[t + 3]), 42, 1);
        require(pl >= 0 && nl >= 0 && pl <= ph && nl <= nh && ph <= ((V)1 << 45) &&
                    nh <= ((V)1 << 45),
                "outward matrix endpoints");
        clo += (__int128)pl - nh;
        chi += (__int128)ph - nl;
        size_t z = (size_t)i * N + j;
        if (!mode) {
          fl[z] = pl;
          fh[z] = ph;
          gl[z] = nl;
          gh[z] = nh;
        } else {
          dpl[z] = pl;
          dnh[z] = nh;
        }
        reset(m);
      }
      if (!mode || i)
        require(clo <= 0 && chi >= 0, "moment conservation");
    }
  }
}
static void first_cells(B M) {
  B eps = xs;
  for (int i = 0; i < N; i++) {
    size_t m = mark();
    B x = xs + i, cut = minimum_b(eps, mul(rh, x)), negpart = neg(S(x, cut)),
      pospart = add(S(x, eps), negpart);
    B er = over(times(mul(mul(M, sq(x)), mul(eps, root(eps))), 2), 3);
    size_t j = (size_t)i * N;
    fl[j] = maxv(0, lattice(sub(pospart, er), 42, 0));
    fh[j] = maxv(0, lattice(add(pospart, er), 42, 1));
    gl[j] = maxv(0, lattice(sub(negpart, er), 42, 0));
    gh[j] = maxv(0, lattice(add(negpart, er), 42, 1));
    if (i) {
      B z = xs + zi[i];
      B dn = neg(sub(S(x, eps), S(z, eps)));
      er = over(times(mul(mul(M, sub(sq(x), sq(z))), mul(eps, root(eps))), 2), 3);
      dpl[j] = 0;
      dnh[j] = maxv(0, lattice(add(dn, er), 42, 1));
    }
    reset(m);
  }
}
static void seeds(B M) {
  size_t m = mark();
  B eps = xs;
  arb_t maxA;
  arb_init(maxA);
  arb_set(maxA, upper(add(aF, over(times(mul(mul(M, eps), mul(k.pi, root(k.pi))), 2), 3))));
  arb_ptr weights = _arb_vec_init(N);
  arb_one(weights);
  for (int i = 1; i < N; i++)
    arb_sqrt(weights + i, xs + i, CAP_PREC);
  for (int i = 1; i < N; i++) {
    size_t q = mark();
    arb_t v;
    arb_init(v);
    for (int j = 0; j < N; j++) {
      size_t z = mark();
      arb_addmul(v, ci(fh[i * N + j]), weights + j, CAP_PREC);
      reset(z);
    }
    arb_div(v, v, xs + i - 1, CAP_PREC);
    arb_max(maxA, maxA, upper(v), CAP_PREC);
    arb_clear(v);
    reset(q);
  }
  show("universal W coefficient", maxA);
  require(lt(maxA, integer(1)), "universal upper seed");
  require(lt(maxA, decimal(".855")),
          "universal W coefficient below 0.855 as stated in Section 8 and Appendix C");
  B d = decimal(".437");
  B seed = sub(over(minimum_b(hs(d, d), hs(integer(1), d)), 2),
               over(mul(mul(M, sub(integer(1), sq(d))), mul(sq(eps), root(eps))), 2));
  show("cusp lower seed", seed);
  require(lt(decimal(".19"), seed), "cusp lower seed");
  for (int i = 0; i < N; i++) {
    size_t q = mark();
    L[i] = maxv(0, lattice(exact_double(raw[(size_t)N * N * 8 + i * 4]), 32, 0));
    U[i] = lattice(minimum_b(root(xs + i), integer(1)), 32, 1);
    reset(q);
  }
  monotone();
  _arb_vec_clear(weights, N);
  arb_clear(maxA);
  reset(m);
}
static B pobj(B p, B hh, B kh) {
  B t = sq(divide(p, hh));
  return add(sub(sq(p), mul(kh, p)),
             add(mul(p, hs(t, delta)), mul(hh, sub(S(integer(1), t), S(delta, t)))));
}
static B dp(B p, B hh, B kh) {
  return add(sub(times(p, 2), kh), hs(sq(divide(p, hh)), delta));
}
static void refine(V *l0out, V *h0out) {
  size_t m = mark();
  V l0 = lattice(decimal(".19"), 32, 0), h0 = BS;
  B eps = xs, err = times(eps, 2);
  arb_ptr fm = _arb_vec_init(N), dm = _arb_vec_init(N), wl = _arb_vec_init(N);
  arb_t kh, bb, fhead;
  arb_init(kh);
  arb_init(bb);
  arb_init(fhead);
  arb_set(kh, minimum_b(hs(delta, delta), hs(integer(1), delta)));
  arb_set(bb, sub(kh, hs(integer(1), delta)));
  arb_set(fhead, sub(S(integer(1), integer(1)), S(integer(1), rh)));
  for (int j = 1; j < N; j++) {
    size_t q = mark();
    B a = divide(xs + j, eps), b = divide(xs + j - 1, eps);
    arb_set(fm + j, sub(S(integer(1), a), S(integer(1), b)));
    arb_set(dm + j, add(sub(fm + j, S(delta, a)), S(delta, b)));
    require(positive(fm + j) && positive(dm + j), "positive cusp tail moments");
    reset(q);
  }
  for (int it = 0; it < 240; it++) {
    size_t q = mark();
    B ll = bi(l0), hh = bi(h0);
    arb_set(wl, ll);
    arb_t hg, outer, ftail;
    arb_init(hg);
    arb_init(outer);
    arb_init(ftail);
    arb_set(hg, hh);
    for (int j = 1; j < N; j++) {
      size_t z = mark();
      V value =
          minv(lattice(wl + j - 1, 32, 0), lattice(divide(bi(L[j - 1]), root(xs + j)), 32, 0));
      arb_set(wl + j, bi(value));
      arb_max(hg, hg, upper(divide(bi(U[j]), root(xs + j - 1))), CAP_PREC);
      arb_addmul(outer, dm + j, wl + j, CAP_PREC);
      arb_addmul(ftail, fm + j, wl + j, CAP_PREC);
      reset(z);
    }
    V nh = minv(h0, lattice(root(add(sub(mul(aF, hg), mul(aG, ll)), err)), 32, 1));
    arb_t pa, pb;
    arb_init(pa);
    arb_init(pb);
    arb_set(pa, lower(mul(ll, root(delta))));
    arb_set(pb, upper(mul(hh, root(delta))));
    if (positive(dp(pa, hh, kh)))
      arb_set(pb, pa);
    else if (negative(dp(pb, hh, kh)))
      arb_set(pa, pb);
    else
      for (int b = 0; b < 65; b++) {
        size_t z = mark();
        B mid = over(add(pa, pb), 2), der = dp(mid, hh, kh);
        if (negative(der))
          arb_set(pa, mid);
        else if (positive(der))
          arb_set(pb, mid);
        else {
          reset(z);
          break;
        }
        reset(z);
      }
    B cc = sub(add(pobj(hull(pa, pb), hh, kh), outer), err);
    V nl = l0;
    if (negative(sub(sub(sq(ll), mul(bb, ll)), cc))) {
      B value = over(add(bb, root(add(sq(bb), times(cc, 4)))), 2);
      nl = maxv(nl, lattice(value, 32, 0));
    }
    B jl = sub(sub(add(mul(fhead, ll), ftail), mul(aG, hh)), err);
    if (positive(jl))
      nl = maxv(nl, lattice(root(jl), 32, 0));
    require(nl <= nh, "cusp refinement ordered");
    l0 = nl;
    h0 = nh;
    L[0] = maxv(L[0], lattice(mul(bi(l0), root(eps)), 32, 0));
    U[0] = minv(U[0], lattice(mul(bi(h0), root(eps)), 32, 1));
    monotone();
    require(advance(N, fl, fh, gl, gh, dpl, dnh, zi, L, U, l0, h0, it < 60 || it % 5 == 0) == 0,
            "exact integer refinement");
    if (it % 10 == 0)
      printf("ITERATION %d cusp_integer=%ld,%ld trough_integer=%ld,%ld scale_bits=32\n", it,
             (long)l0, (long)h0, (long)L[N - 1], (long)U[N - 1]);
    arb_clear(hg);
    arb_clear(outer);
    arb_clear(ftail);
    arb_clear(pa);
    arb_clear(pb);
    reset(q);
  }
  *l0out = l0;
  *h0out = h0;
  show("cusp lower", bi(l0));
  show("cusp upper", bi(h0));
  for (int i = 0; i < N; i++)
    printf("BOUND %d %ld %ld\n", i, (long)L[i], (long)U[i]);
  /* Rounded enclosures printed in Proposition P.uniqueness-enclosures (Section 8). */
  require(lt(decimal(".6105"), bi(l0)), "printed cusp lower bound 0.6105");
  require(lt(bi(h0), decimal(".6399")), "printed cusp upper bound 0.6399");
  require(lt(decimal(".5554"), bi(L[N - 1])), "printed trough lower bound 0.5554");
  require(lt(bi(U[N - 1]), decimal(".6068")), "printed trough upper bound 0.6068");
  _arb_vec_clear(fm, N);
  _arb_vec_clear(dm, N);
  _arb_vec_clear(wl, N);
  arb_clear(kh);
  arb_clear(bb);
  arb_clear(fhead);
  reset(m);
}
static B weight(B x) {
  return add(root(x), times(x, 2));
}
static B fullK(B s) {
  B d = sub(k.pi, ab(sub(ab(s), k.pi)));
  require(positive(d), "comparison avoids pole");
  B val = add(add(divide(k.C, root(d)), divide(zetab(rational(1, 2)), k.pi)),
              add(divide(regular(&k, d, 0), k.pi), divide(integer(1), times(k.pi, 2))));
  B sum = integer(0);
  for (int n = 1; n <= k.N; n++)
    sum = add(sum, mul(k.w + n - 1, cosb(times(d, n))));
  B r = expb(integer(-2)),
    tail = divide(times(power(r, k.N + 1), 2), mul(k.pi, sub(integer(1), r)));
  return inflate(sub(val, divide(sum, k.pi)), tail);
}
static B II(B t) {
  return mul(mul(rational(4, 3), k.C), mul(t, root(t)));
}
static void comparison(V l0, int reject_test) {
  size_t m = mark();
  B eps = xs, M3 = regular_second(&k, times(k.pi, 2), 1);
  show("M3", M3);
  require(lt(M3, decimal(".067")), "M3 bound");
  B Dcoeff = over(times(mul(M3, mul(k.pi, root(k.pi))), 2), 3);
  B t = rh,
    Y = add(sub(add(neg(mul(t, ip(sub(integer(1), t)))), mul(t, ip(add(integer(1), t)))),
                add(II(sub(integer(1), t)), II(add(integer(1), t)))),
            add(neg(times(mul(t, ip(t)), 2)), add(times(II(t), 2), times(II(integer(1)), 2))));
  B absy = sub(mul(rational(8, 3), k.C), times(Y, 2));
  B qsmall = divide(
      add(add(add(k.lambda, times(aG, 2)), times(mul(absy, root(eps)), 2)),
          mul(mul(decimal(".083"), eps), add(over(times(mul(k.pi, root(k.pi)), 2), 3), sq(k.pi)))),
      times(bi(l0), 2));
  show("Q_small", qsmall);
  require(lt(qsmall, decimal(reject_test ? ".5" : ".933")), "cusp comparison target");
  arb_t maxq, nodal;
  arb_init(maxq);
  arb_init(nodal);
  arb_set(maxq, upper(qsmall));
  arb_ptr cs[2] = {_arb_vec_init(13), _arb_vec_init(13)}, weights = _arb_vec_init(N);
  root_coefficients(cs[0], cs[1]);
  arb_set(weights, add(integer(1), times(root(eps), 2)));
  for (int j = 1; j < N; j++) {
    size_t z = mark();
    arb_set(weights + j, weight(xs + j));
    reset(z);
  }
  for (int i = 0; i < N - 1; i++) {
    size_t q = mark();
    B a = xs + i, b = xs + i + 1, cell = hull(a, b);
    arb_ptr low = tmp(), high = tmp(), d = tmp(), c = tmp();
    poly_box(low, d, c, cs[0], 13, lower(a), upper(b));
    poly_box(high, d, c, cs[1], 13, lower(a), upper(b));
    B ratio = hull(lower(low), upper(high));
    require(positive(ratio) && lt(ratio, integer(1)), "root ratio in unit interval");
    B yy = mul(cell, ratio), Sroot = sub(fullK(mul(cell, sub(integer(1), ratio))),
                                         fullK(mul(cell, add(integer(1), ratio))));
    B deriv = add(add(add(k.lambda, mul(Dcoeff, b)), integer(2)),
                  times(mul(absupper(Sroot), absupper(weight(yy))), 2));
    arb_t numerator;
    arb_init(numerator);
    for (int j = 0; j < N; j++) {
      size_t z = mark();
      arb_addmul(numerator, ci(fh[i * N + j] + gh[i * N + j]), weights + j, CAP_PREC);
      reset(z);
    }
    B star = sub(exact_double(raw[(size_t)N * N * 8 + i * 4 + 2]), mul(decimal(".00006"), a)),
      den = mul(weight(a), add(bi(L[i]), star));
    require(positive(den), "comparison denominator");
    B qnode = divide(numerator, den), v = divide(add(numerator, mul(deriv, sub(b, a))), den);
    arb_max(nodal, nodal, upper(qnode), CAP_PREC);
    arb_max(maxq, maxq, upper(v), CAP_PREC);
    printf("COMPARISON_CELL %d ", i);
    arb_printn(v, 22, 0);
    puts("");
    arb_clear(numerator);
    reset(q);
  }
  show("Q_nodal", nodal);
  show("Q_global", maxq);
  require(lt(maxq, decimal(".982")), "full comparison Q<.982");
  require(lt(maxq, decimal(".9815")), "printed comparison bound Q<.9815 (Appendix C)");
  _arb_vec_clear(weights, N);
  _arb_vec_clear(cs[0], 13);
  _arb_vec_clear(cs[1], 13);
  arb_clear(maxq);
  arb_clear(nodal);
  reset(m);
}
int main(int argc, char **argv) {
  require(argc == 3 || argc == 4, "usage: uniqueness grid moments [reject]");
  kernel_init(&k, 72, 36);
  init_homogeneous();
  inputs(argv[1], argv[2]);
  size_t m = mark();
  B M = regular_second(&k, add(k.pi, xs), 0);
  show("M", M);
  require(lt(M, decimal(".083")), "regular second derivative bound");
  require(lt(over(times(mul(M, mul(k.pi, root(k.pi))), 2), 3), integer(2)),
          "cusp error coefficient");
  first_cells(M);
  seeds(M);
  reset(m);
  V l0, h0;
  refine(&l0, &h0);
  comparison(l0, argc == 4);
  printf("PASS uniqueness full cusp and 719 compact cells precision=%d\n", CAP_PREC);
  return 0;
}
