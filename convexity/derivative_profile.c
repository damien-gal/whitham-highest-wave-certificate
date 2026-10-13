/* AI-assisted interval normalized-profile interface, 2026-10-10.
 * Reuses the frozen derivative cusp coefficient and remainder implementation.
 * Compile this translation unit instead of derivative_cusp.c for operator
 * drivers that require both derivative_cusp_errors and this extra API. */
#define DERIVATIVE_CUSP_LIBRARY
#include "derivative_cusp.c"

/* out[j] encloses x^j U^(j)(x)/sqrt(x), j=0,1,2, for all x in X.
 * The source's exact input is fixed at first call, as in derivative_cusp_errors.
 * Returns 1 only for a finite enclosure on 0<X<=0.1 at the working precision CAP_PREC.
 * No mutable global manuscript or existence code is used. */
int derivative_profile_jets(arb_ptr out, const arb_t X,
                            u_params_struct *params, slong prec) {
  arb_t limit, eps, e, a, power, tail, sq;
  arb_init(limit);
  arb_init(eps);
  arb_init(e);
  arb_init(a);
  arb_init(power);
  arb_init(tail);
  arb_init(sq);
  arb_set_str(limit, "0.1", P);
  int ok = prec == P && arb_is_positive(X) && arb_le(X, limit);
  if (!ok) goto done;
  if (!initialized) {
    up = *params;
    build_tables();
  }
  arf_t upper;
  arf_init(upper);
  arb_get_ubound_arf(upper, X, P);
  arb_set_arf(eps, upper);
  arf_clear(upper);
  arb_sqrt(sq, eps, P);
  for (int d = 0; d < 3; ++d) {
    arb_zero(out + d);
    for (int i = 0; i < nu; ++i) {
      exponent(e, us[i].h, us[i].k, us[i].j);
      falling(a, e, d);
      arb_mul(a, a, us[i].c, P);
      exponent(e, us[i].h - 1, us[i].k, us[i].j);
      arb_pow(power, X, e, P);
      arb_addmul(out + d, a, power, P);
    }
    cusp_remainder(tail, eps, d, 0);
    /* Tail/sqrt(x) is a positive increasing majorant. Its endpoint
     * evaluation bounds the whole X cell without division by sqrt(X). */
    arb_div(tail, tail, sq, P);
    arb_add_error(out + d, tail);
    ok = ok && arb_is_finite(out + d);
  }
 done:
  if (!ok) for (int d = 0; d < 3; ++d) arb_indeterminate(out + d);
  arb_clear(limit);
  arb_clear(eps);
  arb_clear(e);
  arb_clear(a);
  arb_clear(power);
  arb_clear(tail);
  arb_clear(sq);
  return ok;
}
