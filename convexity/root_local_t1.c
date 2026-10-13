/* AI-assisted local T1 proof certificate, 2026-10-10.
 * Full mathematical derivation: audit/convexity-wave-2026-10-10/t1-local-proof.md.
 * Shared profile/regular helpers copied from local_t2.c; no inherited T1 norm constants are used. */
#include "constant_functions.h"
#include "derivative_bounds.h"
#include <stdlib.h>
static const slong p=CAP_PREC;
static u_params_struct up;
static void show(const char *s,const arb_t x){printf("%s = ",s);arb_printn(x,25,ARB_STR_MORE);puts("");char *d=arb_dump_str(x);printf("ARB_DUMP %s: %s\n",s,d);flint_free(d);}
/* Independently differentiated complete approximation series and tails. */
static int profile(arb_t B0,arb_t B1,arb_t B2,const arb_t eps){
 arb_ptr B=_arb_vec_init(3); derivative_cusp_errors(B,eps,&up,p);
 arb_set(B0,B);arb_set(B1,B+1);arb_set(B2,B+2);
 int ok=arb_is_nonnegative(B0)&&arb_is_nonnegative(B1)&&arb_is_nonnegative(B2);
 _arb_vec_clear(B,3);return ok;
}
/* Full R(0), and an upper bound for |R''| on |z|<=pi+eps, with
 * K(z)=1/sqrt(2*pi*|z|)+R(z). */
static int regular(arb_t kappa,arb_t M2,const arb_t eps){
 arb_t pi,two,z,r,d,a,b,n,e,tail;arb_init(pi);arb_init(two);arb_init(z);arb_init(r);arb_init(d);arb_init(a);arb_init(b);arb_init(n);arb_init(e);arb_init(tail);
 arb_const_pi(pi,p);arb_mul_ui(two,pi,2,p);arb_add(z,pi,eps,p);arb_div(r,z,two,p);arb_sqr(r,r,p);arb_one(d);arb_sub(d,d,r,p);int ok=arb_is_positive(d);
 /* |zeta(1/2-2m)|/(2m)! <= zeta(5/2)*(2*pi)^(-2m), m>=1.
  * Hence sum of twice differentiated regular Clausen series. */
 arb_set_str(a,"2.5",p);arb_zeta(a,a,p);arb_mul_ui(b,r,3,p);arb_add_ui(b,b,1,p);arb_mul(a,a,b,p);arb_mul_ui(a,a,2,p);arb_pow_ui(b,d,3,p);arb_div(a,a,b,p);arb_sqr(b,two,p);arb_div(a,a,b,p);arb_div(M2,a,pi,p);
 arb_set_str(a,"0.5",p);arb_zeta(kappa,a,p);arb_add(kappa,kappa,a,p);
 for(int i=1;i<=32;i++){
  arb_set_ui(n,i);arb_tanh(a,n,p);arb_sqrt(a,a,p);arb_neg(a,a);arb_add_ui(a,a,1,p);arb_sqrt(b,n,p);arb_div(b,a,b,p);arb_sub(kappa,kappa,b,p);
  arb_mul_ui(b,b,(ulong)i*i,p);arb_div(b,b,pi,p);arb_add(M2,M2,b,p);
 }
 arb_set_str(e,"-0.5",p);compute_tail_error(tail,e,33,p);arb_add_error(kappa,tail);arb_div(kappa,kappa,pi,p);
 arb_set_str(e,"1.5",p);compute_tail_error(tail,e,33,p);arb_div(tail,tail,pi,p);arb_abs(tail,tail);arb_add(M2,M2,tail,p);
 ok=ok&&arb_is_negative(kappa)&&arb_is_positive(M2);
 arb_clear(pi);arb_clear(two);arb_clear(z);arb_clear(r);arb_clear(d);arb_clear(a);arb_clear(b);arb_clear(n);arb_clear(e);arb_clear(tail);return ok;
}
int main(int argc,char**argv){
 populate_parameters(&up,p);
 arb_t eps,pi,lambda,B0,B1,B2,L,E,aE,Amax,cb,kappa,M2,kpi,t,v,bracket,half,r;
 arb_init(eps);arb_init(pi);arb_init(lambda);arb_init(B0);arb_init(B1);arb_init(B2);arb_init(L);arb_init(E);arb_init(aE);arb_init(Amax);arb_init(cb);arb_init(kappa);arb_init(M2);arb_init(kpi);arb_init(t);arb_init(v);arb_init(bracket);arb_init(half);arb_init(r);
 arb_set_str(eps,argc>1?argv[1]:"0.002",p);arb_const_pi(pi,p);arb_div_ui(lambda,pi,8,p);arb_sqrt(lambda,lambda,p);arb_set_str(half,"0.5",p);
 arb_set_str(t,"0.25",p);int ok=arb_equal(t,up.a0k)&&arb_is_positive(eps)&&arb_lt(eps,half)&&arb_gt(up.p0,half);arb_one(t);ok=ok&&arb_lt(up.p0,t);
 if(!ok){puts("LOCAL T1 UNRESOLVED: unsupported cutoff/input");return 2;}
 ok=profile(B0,B1,B2,eps)&&ok;ok=regular(kappa,M2,eps)&&ok;
 arb_pow(E,eps,up.p0,p);arb_mul(t,B0,E,p);arb_sub(L,lambda,t,p);ok=ok&&arb_is_positive(L);
 arb_mul(t,B0,half,p);arb_add(t,t,B1,p);arb_div(aE,t,L,p);
 arb_mul(t,aE,E,p);arb_add(Amax,half,t,p);arb_one(t);ok=ok&&arb_lt(Amax,t);
 Kp(kpi,pi,p);ok=ok&&arb_is_positive(kpi);arb_set_str(cb,"0.997362",p);
 arb_mul(bracket,cb,B0,p);arb_mul_ui(bracket,bracket,2,p);
 arb_mul(t,lambda,aE,p);arb_mul_ui(t,t,16,p);arb_div_ui(t,t,3,p);arb_div(t,t,pi,p);arb_add(bracket,bracket,t,p);
 arb_sub(t,up.p0,half,p);arb_pow(t,eps,t,p);arb_mul(bracket,bracket,t,p);
 arb_sub(t,kappa,kpi,p);arb_neg(v,Amax);arb_add_ui(v,v,2,p);arb_addmul(bracket,t,v,p);
 arb_mul_ui(t,Amax,17,p);arb_div_ui(t,t,12,p);arb_set_ui(v,7);arb_div_ui(v,v,3,p);arb_add(t,t,v,p);arb_mul(t,t,M2,p);arb_sqr(v,eps,p);arb_addmul(bracket,t,v,p);
 show("epsilon",eps);show("B0",B0);show("B1",B1);show("lower normalized u",L);show("aE",aE);show("Amax",Amax);show("kappa",kappa);show("M2",M2);show("K(pi)",kpi);show("signed T1 criterion",bracket);
 arb_set_str(r,"0.00006",p);arb_sqrt(t,eps,p);arb_mul(t,t,r,p);arb_mul_ui(t,t,2,p);arb_add(t,t,bracket,p);show("T1 criterion including radius multiplication",t);
 ok=ok&&arb_is_negative(bracket)&&arb_is_negative(t);
 arb_set_str(r,"-0.196",p);ok=ok&&arb_lt(bracket,r); /* Appendix B prints < -0.196 */
 printf("LOCAL T1 %s: row <= 0.997362 on (0,epsilon]; CB rational premise from existence scalar certificate\n",ok?"PASS":"UNRESOLVED");
 return ok?0:2;
}
