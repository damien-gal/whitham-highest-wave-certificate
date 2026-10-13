#include "constant_functions.h"
static const slong p=CAP_PREC;
static void show(const char *s, const arb_t x) {
  printf("%s = ", s);
  arb_printn(x, 25, ARB_STR_MORE);
  puts("");
  char *dump = arb_dump_str(x);
  printf("ARB_DUMP %s: %s\n", s, dump);
  flint_free(dump);
  fflush(stdout);
}
static void signfunc(arb_t out,const arb_t t){
  arb_t a;
  arb_init(a);
  arb_sub_ui(a,t,1,p);
  arb_rsqrt(out,a,p);
  arb_add_ui(a,t,1,p);
  arb_rsqrt(a,a,p);
  arb_add(out,out,a,p);
  arb_sub_ui(out,out,2,p);
  arb_clear(a);
}
static void primitive(arb_t out,const arb_t t){
  arb_t a,b,c,pi;
  arb_init(a);
  arb_init(b);
  arb_init(c);
  arb_init(pi);
  arb_const_pi(pi,p);
  arb_sub_ui(a,t,1,p);
  arb_sqrt(a,a,p);
  arb_add_ui(b,t,1,p);
  arb_sqrt(b,b,p);
  arb_add(c,a,b,p);
  arb_mul(c,c,t,p);
  arb_mul_ui(c,c,2,p);
  arb_sub(out,a,b,p);
  arb_add(out,out,c,p);
  arb_add_ui(out,out,2,p);
  arb_sqrt(a,t,p);
  arb_mul(a,a,t,p);
  arb_mul(a,a,pi,p);
  arb_mul_ui(a,a,3,p);
  arb_div(out,out,a,p);
  arb_mul_ui(out,out,2,p);
  arb_clear(a);
  arb_clear(b);
  arb_clear(c);
  arb_clear(pi);
}
int main(void){
  u_params_struct up;
  populate_parameters(&up,p);
  arb_t a,b,m,t,cb,eps,pi,A,D,q,r,beta,invariance,contraction,nu,mu,lo;
  arb_init(a);
  arb_init(b);
  arb_init(m);
  arb_init(t);
  arb_init(cb);
  arb_init(eps);
  arb_init(pi);
  arb_init(A);
  arb_init(D);
  arb_init(q);
  arb_init(r);
  arb_init(beta);
  arb_init(invariance);
  arb_init(contraction);
  arb_init(nu);
  arb_init(mu);
  arb_init(lo);
  arb_const_pi(pi,p);
  arb_set_str(a,"1.5",p);
  arb_set_ui(b,2);
  signfunc(t,a);
  int ok=arb_is_positive(t);
  signfunc(t,b);
  ok=ok&&arb_is_negative(t);
  for(int k=0;k<90;k++){
    arb_add(m,a,b,p);
    arb_mul_2exp_si(m,m,-1);
    signfunc(t,m);
    if(arb_is_positive(t))arb_set(a,m);
    else if(arb_is_negative(t))arb_set(b,m);
    else{
      ok=0;
      break;
    }
  }
  arb_union(m,a,b,p);
  primitive(cb,m);
  arb_mul_ui(cb,cb,2,p);
  arb_mul_ui(t,pi,3,p);
  arb_ui_div(t,8,t,p);
  arb_sub(cb,cb,t,p);
  show("C_B",cb);
  arb_one(t);
  ok=ok&&arb_lt(cb,t);
  arb_set_str(eps,"0.05",p);
  compute_c_hat_minus(t,eps,7,&up,p);
  show("c_hat_minus",t);
  ok=ok&&arb_is_positive(t);
  compute_c_hat_plus(a,eps,7,&up,p);
  show("c_hat_plus",a);
  ok=ok&&arb_is_positive(a);
  arb_pow(t,eps,up.p0,p);
  arb_mul(t,t,a,p);
  arb_div_ui(lo,pi,8,p);
  arb_sqrt(lo,lo,p);
  arb_sub(lo,lo,t,p);
  show("hat denominator lower bound",lo);
  ok=ok&&arb_is_positive(lo);
  compute_validation_T0(t,eps,7,&up,p);
  show("local T0 scalar residual",t);
  ok=ok&&arb_is_negative(t);
  compute_cT0_1(t,eps,&up,p);
  ok=ok&&arb_is_positive(t);
  compute_cT0_2(t,eps,&up,p);
  ok=ok&&arb_is_positive(t);
  arb_set_str(t,"0.5",p);
  ok=ok&&arb_gt(up.p0,t);
  arb_one(t);
  ok=ok&&arb_lt(up.p0,t);
  arb_set_str(A,"2.7",p);
  arb_set_str(D,"0.0000001",p);
  arb_set_str(q,"0.998",p);
  ok=ok&&arb_lt(cb,q);
  arb_set_str(r,"0.00006",p);
  arb_neg(beta,q);
  arb_add_ui(beta,beta,1,p);
  arb_inv(beta,beta,p);
  /* Local alpha <= sqrt(eps)/(2 lower(u/sqrt x)). */
  arb_sqrt(t,eps,p);
  arb_div(t,t,lo,p);
  arb_mul_2exp_si(t,t,-1);
  show("cusp alpha upper bound",t);
  ok=ok&&arb_lt(t,A);
  arb_sqr(t,r,p);
  arb_mul(t,t,A,p);
  arb_add(t,t,D,p);
  arb_mul(t,t,beta,p);
  arb_sub(invariance,r,t,p);
  arb_mul(contraction,beta,A,p);
  arb_mul(contraction,contraction,r,p);
  arb_mul_ui(contraction,contraction,2,p);
  show("A",A);
  show("delta",D);
  show("q",q);
  show("radius",r);
  show("invariance margin",invariance);
  show("contraction factor",contraction);
  arb_one(t);
  ok=ok&&arb_is_positive(invariance)&&arb_lt(contraction,t);
  arb_mul(t,A,r,p);
  arb_mul_ui(t,t,2,p);
  arb_neg(t,t);
  arb_add_ui(t,t,1,p);
  show("positivity factor 1-2*A*r",t);
  ok=ok&&arb_is_positive(t);
  /* Inherited nu already includes radius r_old=0.0000359. Enlarge the
    symmetric enclosure by 4*pi*(r-r_old); no double counting. */
  compute_nu(nu,&up,p);
  arb_set_str(t,"0.0000359",p);
  arb_sub(t,r,t,p);
  ok=ok&&arb_is_nonnegative(t);
  arb_mul(t,t,pi,p);
  arb_mul_ui(t,t,4,p);
  arb_add_error(nu,t);
  show("nu = 4 L u(0)",nu);
  arb_neg(t,nu);
  arb_add_ui(t,t,1,p);
  show("speed radicand",t);
  ok=ok&&arb_is_positive(t)&&arb_is_positive(nu);
  arb_sqrt(mu,t,p);
  arb_neg(mu,mu);
  arb_add_ui(mu,mu,1,p);
  show("lower-branch speed",mu);
  /* Printed enclosures (Theorem 1.1 and Proposition P.mu): the recorded
     balls for nu and mu must lie inside the manuscript's bounds. */
  arb_set_str(t,"0.9455",p);ok=ok&&arb_gt(nu,t);
  arb_set_str(t,"0.9471",p);ok=ok&&arb_lt(nu,t);
  arb_set_str(t,"0.766",p);ok=ok&&arb_gt(mu,t);
  arb_set_str(t,"0.770",p);ok=ok&&arb_lt(mu,t);
  /* C_B premise hard-coded by the local T1/T2 stages and Appendix A. */
  arb_set_str(t,"0.997362",p);ok=ok&&arb_lt(cb,t);
  printf("SCALAR CONDITIONS %s (conditional on matching full-domain certificates)\n",ok?"PASS":"UNRESOLVED");
  return ok?0:2;
}
