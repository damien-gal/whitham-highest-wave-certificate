/* kernel.h: ball arithmetic for the Whitham kernel K, used by the kernel-sign
 * certificate (sign.c), the regression check (test_kernel.c) and the
 * uniqueness comparison (uniqueness.c). See kernel.c for the formulas.
 */
#ifndef WHITHAM_KERNEL_H
#define WHITHAM_KERNEL_H
#include "ball.h"
typedef struct {
  int M, N;
  arb_t pi, C, lambda, A;
  arb_ptr coeff, w;
} Kernel;
void kernel_init(Kernel *k, int M, int N);
void kernel_clear(Kernel *k);
arb_ptr regular(Kernel *k, B z, int derivative);
arb_ptr homogeneous(Kernel *k, B t);
arb_ptr curve(Kernel *k, B x, B t);
arb_ptr curve_derivative(Kernel *k, B x, B t, B tp);
arb_ptr small_curve(Kernel *k, B t, B eps);
arb_ptr fourier_second(Kernel *k);
arb_ptr regular_second(Kernel *k, B xmax, int images);
void root_coefficients(arb_ptr minus, arb_ptr plus);
void poly_box(arb_t value, arb_t deriv, arb_t center, arb_srcptr cs, int n, B a, B b);
#endif
