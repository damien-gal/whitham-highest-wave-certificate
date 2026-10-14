/* closure.c: scalar assembly of the correction bounds (stage closure).
 *
 * Combines the exact targets of inputs/closure.json (contraction factors,
 * defect bounds, coefficient bounds and curvature budget) into the radii of
 * the fixed-point argument of Sections 6-7, with ball arithmetic at CAP_PREC.
 * The uniform bounds it consumes are premises that the full-domain stages
 * certify separately; run.py supplies their conjunction.
 */
#include "ball.h"
int main(void) {
  B A = decimal("2.7"), r = decimal(".00006"), q0 = decimal(".998"), d0 = decimal(".0000001"),
    d1 = decimal(".0000005"), d2 = decimal(".00001"), a1 = decimal(".32"), a2 = decimal("1.383"),
    af = decimal(".414");
  B beta0 = divide(integer(1), sub(integer(1), q0));
  B shift = times(mul(A, r), 2);
  B den = sub(sub(integer(1), q0), shift);
  require(positive(den), "derivative contraction margin");
  B e1 = divide(add(d1, mul(a1, sq(r))), den);
  B e2 = divide(add(add(add(add(d2, mul(af, e1)), mul(a2, sq(r))), times(mul(A, sq(e1)), 2)),
                    mul(mul(times(a1, 4), r), e1)),
                den);
  B invariance = sub(r, mul(beta0, add(d0, mul(A, sq(r))))),
    contraction = times(mul(mul(beta0, A), r), 2), curvature = sub(decimal(".089"), e2);
  require(positive(invariance) && lt(contraction, integer(1)) && positive(sub(integer(1), shift)),
          "existence closure");
  require(lt(e1, decimal(".0003")) && lt(e2, decimal(".08014")) && lt(decimal(".00886"), curvature),
          "derivative/convexity closure");
  show("epsilon1", e1);
  show("epsilon2", e2);
  show("curvature margin", curvature);
  show("invariance margin", invariance);
  show("existence contraction", contraction);
  show("derivative contraction", add(q0, shift));
  puts("PASS scalar closure CONDITIONAL on matching full-domain covers");
  return 0;
}
