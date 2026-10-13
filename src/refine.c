/* refine.c: exact integer refinement of the a priori enclosures (Section 8 and
 * Appendix C of the paper).
 *
 * Bounds are stored at scale 2^32 and the nonnegative matrix entries at scale
 * 2^42. Products and sums use signed 128-bit integers; the inputs are capped
 * below 2^45 for coefficients and 2^33 for bounds, so that N <= 10000 keeps
 * every sum below 2^93. No floating-point arithmetic occurs in this file.
 * test_refine.c checks the discrete minimization exhaustively on small domains.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 I;
typedef int64_t V;
#define SHIFT 10
static I square(V q) {
  return ((I)q * q) << SHIFT;
}
typedef struct {
  int n;
  const V *v;
  V *s;
  I *sv;
} Sum;
static Sum make(const V *a, const V *v, int n) {
  Sum r = {n, v, malloc((n + 1) * sizeof(V)), malloc((n + 1) * sizeof(I))};
  r.s[0] = 0;
  r.sv[0] = 0;
  for (int j = 0; j < n; j++) {
    r.s[j + 1] = r.s[j] + a[j];
    r.sv[j + 1] = r.sv[j] + (I)a[j] * v[j];
  }
  return r;
}
static void clear(Sum *s) {
  free(s->s);
  free(s->sv);
}
static int split(const Sum *s, V p) {
  int a = 0, b = s->n;
  while (a < b) {
    int m = (a + b) / 2;
    if (s->v[m] < p)
      a = m + 1;
    else
      b = m;
  }
  return a;
}
static I minimum(const Sum *s, V p) {
  int j = split(s, p);
  return s->sv[j] + (I)p * (s->s[s->n] - s->s[j]);
}
static I maximum(const Sum *s, V p) {
  int j = split(s, p);
  return (I)p * s->s[j] + s->sv[s->n] - s->sv[j];
}
static V isqrt_up(I v, V top) {
  V lo = 0, hi = top;
  while (lo < hi) {
    V m = lo + (hi - lo) / 2;
    if (square(m) >= v)
      hi = m;
    else
      lo = m + 1;
  }
  return lo;
}
static V isqrt_down(I v, V top) {
  V q = isqrt_up(v, top);
  if (square(q) > v && q)
    q--;
  return q;
}
typedef struct {
  Sum a, b, m, c;
  V pu, tu;
  I err, base, outer;
  int difference;
} Row;
static I pobj(const Row *r, V p) {
  return square(p) - minimum(&r->b, p);
}
static I inner(const Row *r, V t) {
  if (!r->difference)
    return maximum(&r->a, t) - minimum(&r->b, t);
  V p = t < r->pu ? t : r->pu;
  return pobj(r, p) + maximum(&r->a, t) - minimum(&r->m, t);
}
static V minp(const Row *r, V lo, V hi) {
  while (lo < hi) {
    V m = lo + (hi - lo) / 2;
    if (pobj(r, m + 1) < pobj(r, m))
      lo = m + 1;
    else
      hi = m;
  }
  return lo;
}
static V mint(const Row *r, V lo, V hi) {
  while (lo < hi) {
    V m = lo + (hi - lo) / 2;
    if (inner(r, m + 1) < inner(r, m))
      lo = m + 1;
    else
      hi = m;
  }
  return lo;
}
static I residual(const Row *r, V q) {
  V t = q < r->tu ? q : r->tu;
  return square(q) - inner(r, t) + r->err - r->base - r->outer;
}
static V improve(const Row *r, V old, V upper) {
  if (residual(r, old) >= 0)
    return old;
  /* Largest integer b with strict whole-interval exclusion. */
  V a = old, b = upper;
  while (a < b) {
    V m = a + (b - a + 1) / 2;
    if (residual(r, m) < 0)
      a = m;
    else
      b = m - 1;
  }
  return a;
}
int advance(int n, const V *fl, const V *fh, const V *gl, const V *gh, const V *dpl, const V *dnh,
            const int *zi, V *L, V *U, V l0, V h0, int differences) {
  if (n < 2 || n > 10000 || l0 < 0 || h0 < l0 || h0 > ((V)1 << 33))
    return -1;
  for (int i = 0; i < n * n; i++)
    if (fl[i] < 0 || gl[i] < 0 || fh[i] < fl[i] || gh[i] < gl[i] || fh[i] > ((V)1 << 45) ||
        gh[i] > ((V)1 << 45) || dpl[i] < 0 || dnh[i] < 0 || dpl[i] > ((V)1 << 45) ||
        dnh[i] > ((V)1 << 45))
      return -2;
  for (int i = 0; i < n; i++)
    if (L[i] < 0 || U[i] < L[i] || U[i] > ((V)1 << 33) ||
        (i && (L[i] < L[i - 1] || U[i] < U[i - 1])))
      return -3;
  V *lc = malloc(n * sizeof(V)), *uc = malloc(n * sizeof(V)), *nl = malloc(n * sizeof(V)),
    *nu = malloc(n * sizeof(V));
  lc[0] = l0;
  uc[0] = h0;
  for (int j = 1; j < n; j++) {
    lc[j] = L[j - 1];
    uc[j] = U[j];
  }
  for (int i = 0; i < n; i++) {
    I low = 0, high = 0;
    for (int j = 0; j < n; j++) {
      low += (I)fl[i * n + j] * lc[j] - (I)gh[i * n + j] * uc[j];
      high += (I)fh[i * n + j] * uc[j] - (I)gl[i * n + j] * lc[j];
    }
    if (high < 0)
      return -4;
    nl[i] = L[i];
    nu[i] = U[i];
    if (low > 0) {
      V q = isqrt_down(low, U[i]);
      if (q > nl[i])
        nl[i] = q;
    }
    if (high < square(U[i]))
      nu[i] = isqrt_up(high, U[i]);
    if (nu[i] < nl[i])
      return -5;
  }
  for (int i = 1; i < n; i++) {
    const V *fp = fl + i * n, *gn = gh + i * n;
    Row r = {0};
    r.a = make(fp + 1, lc + 1, i);
    r.b = make(gn + 1, uc + 1, i);
    r.c = make(fp + i + 1, lc + i + 1, n - i - 1);
    r.base = (I)fp[0] * l0 - (I)gn[0] * h0;
    r.outer = maximum(&r.c, L[i]);
    r.err = r.a.s[r.a.n] + r.b.s[r.b.n];
    r.tu = mint(&r, 0, U[i]);
    V q = improve(&r, nl[i], nu[i]);
    if (q > nl[i])
      nl[i] = q;
    clear(&r.a);
    clear(&r.b);
    clear(&r.c);
    if (differences) {
      int k = zi[i];
      if (k < 0 || k >= i)
        return -6;
      fp = dpl + i * n;
      gn = dnh + i * n;
      memset(&r, 0, sizeof(r));
      r.difference = 1;
      r.a = make(fp + k + 1, lc + k + 1, i - k);
      r.b = make(gn + 1, uc + 1, k);
      r.m = make(gn + k + 1, uc + k + 1, i - k);
      r.c = make(fp + i + 1, lc + i + 1, n - i - 1);
      r.base = (I)fp[0] * l0 - (I)gn[0] * h0;
      r.outer = maximum(&r.c, L[i]);
      r.pu = minp(&r, L[k], U[k]);
      r.tu = mint(&r, L[k], U[i]);
      r.err = ((I)2 * U[k] << SHIFT) + r.b.s[r.b.n] + r.a.s[r.a.n] + r.m.s[r.m.n];
      q = improve(&r, nl[i], nu[i]);
      if (q > nl[i])
        nl[i] = q;
      clear(&r.a);
      clear(&r.b);
      clear(&r.m);
      clear(&r.c);
    }
  }
  for (int i = 1; i < n; i++)
    if (nl[i] < nl[i - 1])
      nl[i] = nl[i - 1];
  for (int i = n - 2; i >= 0; i--)
    if (nu[i] > nu[i + 1])
      nu[i] = nu[i + 1];
  for (int i = 0; i < n; i++)
    if (nl[i] > nu[i])
      return -7;
  memcpy(L, nl, n * sizeof(V));
  memcpy(U, nu, n * sizeof(V));
  free(lc);
  free(uc);
  free(nl);
  free(nu);
  return 0;
}
