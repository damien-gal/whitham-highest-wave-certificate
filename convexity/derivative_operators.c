/* AI side-wave verifier. Full kernel, zero mode and both periodic poles. */
#include "basic_functions.h"
#include "complex_functions.h"
#include "acb_calc.h"
#include <stdlib.h>
#define M 24
static int analytic=0;
static slong p=CAP_PREC;
static arb_t pi,L,c,rootL,k[M],kl[M],bin[M][M],s[33],p10,p20,tp,lim;
static u_params_struct up;
static long calls,nodes,leaves;
int derivative_profile_jets(arb_ptr,const arb_t,u_params_struct*,slong);
typedef struct{
  arb_t x,u,a,d,f,w,A,B,H;
  int region,op;
}
CT;
#define AI(x) arb_init(x)
#define CI(x) acb_init(x)
#define AC(x) arb_clear(x)
#define CC(x) acb_clear(x)
static void setup(void){
  arb_t t,z,b,h,e;
  AI(t);
  AI(z);
  AI(b);
  AI(h);
  AI(e);
  AI(pi);
  AI(L);
  AI(c);
  AI(rootL);
  AI(p10);
  AI(p20);
  AI(tp);
  AI(lim);
  arb_const_pi(pi,p);
  arb_mul_ui(L,pi,2,p);
  arb_sqrt(rootL,L,p);
  arb_inv(c,rootL,p);
  for(int j=0;j<M;j++){
    AI(k[j]);
    AI(kl[j]);
  }
  arb_set_str(t,".5",p);
  arb_zeta(t,t,p);
  arb_sqrt_ui(z,2,p);
  arb_sub_ui(z,z,1,p);
  arb_mul(t,t,z,p);
  arb_set_str(z,".5",p);
  arb_add(t,t,z,p);
  arb_div(t,t,pi,p);
  arb_sqrt(z,pi,p);
  arb_div(z,c,z,p);
  arb_mul_ui(z,z,2,p);
  arb_sub(k[0],t,z,p);
  arb_set_str(t,".5",p);
  arb_zeta(t,t,p);
  arb_set_str(z,".5",p);
  arb_add(t,t,z,p);
  arb_div(kl[0],t,pi,p);
  arb_one(b);
  for(int j=1;j<M;j++){
    arb_mul_ui(b,b,4*j-3,p);
    arb_div_ui(b,b,4*j-2,p);
    arb_mul_ui(b,b,4*j-1,p);
    arb_div_ui(b,b,4*j,p);
    arb_set_ui(e,2*j);
    arb_set_str(z,".5",p);
    arb_add(e,e,z,p);
    arb_set_str(h,"1.5",p);
    arb_hurwitz_zeta(t,e,h,p);
    arb_pow(z,L,e,p);
    arb_div(t,t,z,p);
    arb_mul(t,t,b,p);
    arb_mul(t,t,c,p);
    arb_mul_ui(k[j],t,2,p);
    arb_zeta(t,e,p);
    arb_div(t,t,z,p);
    arb_mul(t,t,b,p);
    arb_mul(t,t,c,p);
    arb_mul_ui(kl[j],t,2,p);
  }
  for(int j=0;j<M;j++)for(int m=0;m<=j;m++){
    AI(bin[j][m]);
    arb_one(t);
    for(int q=1;q<=2*m;q++){
      arb_mul_ui(t,t,2*j-q+1,p);
      arb_div_ui(t,t,q,p);
    }
    arb_mul(bin[j][m],t,k[j],p);
  }
  arb_zero(p10);
  arb_zero(p20);
  for(int j=0;j<M;j++){
    arb_pow_ui(t,pi,2*j+1,p);
    arb_mul(t,t,k[j],p);
    arb_div_ui(t,t,2*j+1,p);
    arb_sub(p10,p10,t,p);
    arb_mul(t,t,pi,p);
    arb_div_ui(t,t,2*j+2,p);
    arb_add(p20,p20,t,p);
  }
  for(int n=1;n<=32;n++){
    AI(s[n]);
    arb_set_ui(t,n);
    arb_tanh(t,t,p);
    arb_sqrt(t,t,p);
    arb_sub_ui(t,t,1,p);
    arb_sqrt_ui(z,n,p);
    arb_mul(z,z,pi,p);
    arb_div(s[n],t,z,p);
  }
  arb_mul_ui(t,pi,3,p);
  arb_sqrt(t,t,p);
  arb_div(tp,c,t,p);
  arb_mul_ui(tp,tp,2,p);
  arb_set_ui(t,4*M-1);
  arb_ui_div(t,3,t,p);
  arb_add_ui(t,t,1,p);
  arb_mul(tp,tp,t,p);
  populate_parameters(&up,p);
  AC(t);
  AC(z);
  AC(b);
  AC(h);
  AC(e);
}
static void tail(arb_t e,const arb_t rad){
  arb_t r,t;
  AI(r);
  AI(t);
  arb_mul_ui(r,pi,3,p);
  arb_div(r,rad,r,p);
  arb_sqr(r,r,p);
  arb_one(t);
  arb_sub(t,t,r,p);
  if(!arb_is_positive(t))arb_indeterminate(e);
  else{
    arb_pow_ui(e,r,M,p);
    arb_mul(e,e,tp,p);
    arb_div(e,e,t,p);
  }
  AC(r);
  AC(t);
}
static void local_kernel(acb_t K,acb_t I,acb_t J,const acb_t z,const acb_ptr f0,int omit){
  acb_t z2,rt,t,v,si,co;
  arb_t a,r,e,d;
  CI(z2);
  CI(rt);
  CI(t);
  CI(v);
  CI(si);
  CI(co);
  AI(a);
  AI(r);
  AI(e);
  AI(d);
  acb_sqr(z2,z,p);
  if(f0)acb_set(rt,f0);
  else acb_sqrt_analytic(rt,z,analytic,p);
  acb_zero(K);
  acb_zero(I);
  acb_zero(J);
  for(int j=M-1;j>=0;j--){
    acb_mul(K,K,z2,p);
    acb_add_arb(K,K,kl[j],p);
    arb_div_ui(a,kl[j],2*j+1,p);
    acb_mul(I,I,z2,p);
    acb_add_arb(I,I,a,p);
    arb_div_ui(a,a,2*j+2,p);
    acb_mul(J,J,z2,p);
    acb_add_arb(J,J,a,p);
  }
  acb_mul(I,I,z,p);
  acb_mul(J,J,z2,p);
  if(!(omit&1)){
    acb_inv(v,rt,p);
    acb_mul_arb(v,v,c,p);
    acb_add(K,K,v,p);
  }
  acb_mul_arb(v,rt,c,p);
  acb_mul_ui(v,v,2,p);
  acb_add(I,I,v,p);
  acb_mul(v,v,z,p);
  acb_mul_ui(v,v,2,p);
  acb_div_ui(v,v,3,p);
  acb_add(J,J,v,p);
  for(int n=1;n<=32;n++){
    acb_mul_ui(t,z,n,p);
    acb_sin_cos(si,co,t,p);
    acb_mul_arb(v,co,s[n],p);
    acb_add(K,K,v,p);
    acb_mul_arb(v,si,s[n],p);
    acb_div_ui(v,v,n,p);
    acb_add(I,I,v,p);
    acb_mul_2exp_si(t,t,-1);
    acb_sin(t,t,p);
    acb_sqr(t,t,p);
    acb_mul_arb(t,t,s[n],p);
    acb_mul_ui(t,t,2,p);
    acb_div_ui(t,t,n*n,p);
    acb_add(J,J,t,p);
  }
  acb_abs(a,z,p);
  arb_div(r,a,L,p);
  arb_sqr(r,r,p);
  arb_one(d);
  arb_sub(d,d,r,p);
  arb_pow_ui(e,r,M,p);
  arb_mul(e,e,c,p);
  arb_mul_ui(e,e,4,p);
  arb_div(e,e,rootL,p);
  arb_div(e,e,d,p);
  arb_abs(r,acb_imagref(z));
  arb_set_ui(d,2);
  arb_sub(d,d,r,p);
  if(!arb_is_positive(d))arb_indeterminate(e);
  else{
    arb_neg(r,d);
    arb_exp(r,r,p);
    arb_one(d);
    arb_sub(d,d,r,p);
    arb_pow_ui(r,r,33,p);
    arb_mul_ui(r,r,4,p);
    arb_div(r,r,d,p);
    arb_div(r,r,pi,p);
  }
  arb_add(d,e,r,p);
  acb_add_error_arb(K,d);
  arb_mul(e,e,a,p);
  arb_add(d,e,r,p);
  acb_add_error_arb(I,d);
  arb_mul(e,e,a,p);
  arb_add(d,e,r,p);
  acb_add_error_arb(J,d);
  /* When caller omits only the far pole, this representation must subtract it. */  if(omit&2){
    acb_neg(t,z);
    acb_add_arb(t,t,L,p);
    acb_rsqrt_analytic(t,t,analytic,p);
    acb_mul_arb(t,t,c,p);
    acb_sub(K,K,t,p);
  }
  CC(z2);
  CC(rt);
  CC(t);
  CC(v);
  CC(si);
  CC(co);
  AC(a);
  AC(r);
  AC(e);
  AC(d);
}
/* z in [0,2pi], optional forced roots for desingularization. */ static void kernel(acb_t K,acb_t I,acb_t J,const acb_t z,const acb_ptr f0,const acb_ptr fL,int omit){
  arb_t test,bound;
  AI(test);
  AI(bound);
  acb_abs(test,z,p);
  arb_set_ui(bound,4);
  int local=arb_lt(test,bound);
  AC(test);
  AC(bound);
  if(local){
    local_kernel(K,I,J,z,f0,omit);
    return;
  }
  acb_t w,v2,r0,rL,t,v,P,P1,P2,co,si;
  arb_t a,e,e0,r,ie,je,ex;
  CI(w);
  CI(v2);
  CI(r0);
  CI(rL);
  CI(t);
  CI(v);
  CI(P);
  CI(P1);
  CI(P2);
  CI(co);
  CI(si);
  AI(a);
  AI(e);
  AI(e0);
  AI(r);
  AI(ie);
  AI(je);
  AI(ex);
  acb_sub_arb(w,z,pi,p);
  acb_sqr(v2,w,p);
  if(f0)acb_set(r0,f0);
  else acb_sqrt_analytic(r0,z,analytic,p);
  acb_neg(t,z);
  acb_add_arb(t,t,L,p);
  if(fL)acb_set(rL,fL);
  else acb_sqrt_analytic(rL,t,analytic,p);
  acb_zero(P);
  acb_zero(P1);
  acb_zero(P2);
  for(int j=M-1;j>=0;j--){
    acb_mul(P,P,v2,p);
    acb_add_arb(P,P,k[j],p);
    arb_div_ui(a,k[j],2*j+1,p);
    acb_mul(P1,P1,v2,p);
    acb_add_arb(P1,P1,a,p);
    arb_div_ui(a,a,2*j+2,p);
    acb_mul(P2,P2,v2,p);
    acb_add_arb(P2,P2,a,p);
  }
  acb_mul(P1,P1,w,p);
  acb_mul(P2,P2,v2,p);
  acb_set(K,P);
  if(!(omit&1)){
    acb_inv(v,r0,p);
    acb_mul_arb(v,v,c,p);
    acb_add(K,K,v,p);
  }
  if(!(omit&2)){
    acb_inv(v,rL,p);
    acb_mul_arb(v,v,c,p);
    acb_add(K,K,v,p);
  }
  acb_sub(I,r0,rL,p);
  acb_add_arb(I,I,rootL,p);
  acb_mul_arb(I,I,c,p);
  acb_mul_ui(I,I,2,p);
  acb_add(I,I,P1,p);
  acb_sub_arb(I,I,p10,p);
  acb_mul(J,z,r0,p);
  acb_mul(t,t,rL,p);
  acb_add(J,J,t,p);
  arb_mul(a,L,rootL,p);
  acb_sub_arb(J,J,a,p);
  acb_mul_arb(J,J,c,p);
  acb_mul_ui(J,J,4,p);
  acb_div_ui(J,J,3,p);
  arb_mul(a,c,rootL,p);
  arb_mul_ui(a,a,2,p);
  acb_mul_arb(v,z,a,p);
  acb_add(J,J,v,p);
  acb_add(J,J,P2,p);
  acb_sub_arb(J,J,p20,p);
  acb_mul_arb(v,z,p10,p);
  acb_sub(J,J,v,p);
  for(int n=1;n<=32;n++){
    acb_mul_ui(t,z,n,p);
    acb_sin_cos(si,co,t,p);
    acb_mul_arb(v,co,s[n],p);
    acb_add(K,K,v,p);
    acb_mul_arb(v,si,s[n],p);
    acb_div_ui(v,v,n,p);
    acb_add(I,I,v,p);
    acb_one(t);
    acb_sub(t,t,co,p);
    acb_mul_arb(v,t,s[n],p);
    acb_div_ui(v,v,n*n,p);
    acb_add(J,J,v,p);
  }
  acb_abs(r,w,p);
  tail(e,r);
  tail(e0,pi);
  arb_mul(ie,e,r,p);
  arb_mul(a,e0,pi,p);
  arb_add(ie,ie,a,p);
  arb_sqr(a,r,p);
  arb_mul(je,e,a,p);
  arb_sqr(a,pi,p);
  arb_mul(a,a,e0,p);
  arb_add(je,je,a,p);
  acb_abs(a,z,p);
  arb_mul(a,a,pi,p);
  arb_mul(a,a,e0,p);
  arb_add(je,je,a,p);
  arb_abs(ex,acb_imagref(z));
  arb_set_ui(a,2);
  arb_sub(a,a,ex,p);
  if(!arb_is_positive(a))arb_indeterminate(e);
  else{
    arb_neg(ex,a);
    arb_exp(ex,ex,p);
    arb_one(a);
    arb_sub(a,a,ex,p);
    arb_pow_ui(ex,ex,33,p);
    arb_mul_ui(ex,ex,4,p);
    arb_div(ex,ex,a,p);
    arb_div(ex,ex,pi,p);
    arb_add(e,e,ex,p);
    arb_add(ie,ie,ex,p);
    arb_add(je,je,ex,p);
  }
  acb_add_error_arb(K,e);
  acb_add_error_arb(I,ie);
  acb_add_error_arb(J,je);
  CC(w);
  CC(v2);
  CC(r0);
  CC(rL);
  CC(t);
  CC(v);
  CC(P);
  CC(P1);
  CC(P2);
  CC(co);
  CC(si);
  AC(a);
  AC(e);
  AC(e0);
  AC(r);
  AC(ie);
  AC(je);
  AC(ex);
}
/* Stable outer raw kernel from its even Taylor series in x. */ static void outer_raw(acb_t out,const acb_t y,CT*o){
  acb_t z,r,r2,pow,weight,tmp,v,yy,reg;
  arb_t A,B,b,co,tx,nx,sn,cn,e,rad,aa,bb,rr,fac;
  CI(z);
  CI(r);
  CI(r2);
  CI(pow);
  CI(weight);
  CI(tmp);
  CI(v);
  CI(yy);
  CI(reg);
  AI(A);
  AI(B);
  AI(b);
  AI(co);
  AI(tx);
  AI(nx);
  AI(sn);
  AI(cn);
  AI(e);
  AI(rad);
  AI(aa);
  AI(bb);
  AI(rr);
  AI(fac);
  arb_set(A,o->A);
  arb_set(B,o->B);
  acb_zero(out);
  for(int pole=0;pole<2;pole++){
    if(!pole)acb_set(z,y);
    else{
      acb_neg(z,y);
      acb_add_arb(z,z,L,p);
    }
    acb_set_arb(r,o->x);
    acb_div(r,r,z,p);
    acb_sqr(r2,r,p);
    acb_one(pow);
    acb_zero(v);
    arb_one(b);
    for(int m=0;m<M;m++){
      arb_div_ui(co,A,2*m+1,p);
      arb_mul_ui(co,co,2,p);
      arb_neg(co,co);
      arb_add_ui(co,co,1,p);
      arb_div_ui(tx,B,(2*m+1)*(2*m+2),p);
      arb_add(co,co,tx,p);
      arb_mul(co,co,b,p);
      acb_mul_arb(tmp,pow,co,p);
      acb_add(v,v,tmp,p);
      acb_mul(pow,pow,r2,p);
      arb_mul_ui(b,b,4*m+1,p);
      arb_div_ui(b,b,4*m+2,p);
      arb_mul_ui(b,b,4*m+3,p);
      arb_div_ui(b,b,4*m+4,p);
    }
    acb_rsqrt_analytic(tmp,z,analytic,p);
    acb_mul_arb(tmp,tmp,c,p);
    acb_mul_ui(tmp,tmp,2,p);
    acb_mul(v,v,tmp,p);
    acb_add(out,out,v,p);
    acb_abs(rr,r2,p);
    arb_one(aa);
    arb_sub(aa,aa,rr,p);
    if(!arb_is_positive(aa)){
      arb_indeterminate(e);
    }
    else {
      arb_pow_ui(e,rr,M,p);
      arb_div(e,e,aa,p);
    }
    arb_abs(aa,A);
    arb_mul_ui(aa,aa,2,p);
    arb_abs(bb,B);
    arb_add(aa,aa,bb,p);
    arb_add_ui(aa,aa,1,p);
    arb_mul(e,e,aa,p);
    acb_abs(aa,tmp,p);
    arb_mul(e,e,aa,p);
    if(!arb_is_positive(aa)||!arb_is_finite(e))arb_indeterminate(e);
    acb_add_error_arb(out,e);
  }
  acb_sub_arb(yy,y,pi,p);
  acb_sqr(yy,yy,p);
  arb_one(tx);
  acb_zero(reg);
  for(int m=0;m<M;m++){
    acb_zero(v);
    for(int j=M-1;j>=m;j--){
      acb_mul(v,v,yy,p);
      acb_add_arb(v,v,bin[j][m],p);
    }
    arb_div_ui(co,A,2*m+1,p);
    arb_mul_ui(co,co,2,p);
    arb_neg(co,co);
    arb_add_ui(co,co,1,p);
    arb_div_ui(b,B,(2*m+1)*(2*m+2),p);
    arb_add(co,co,b,p);
    arb_mul(co,co,tx,p);
    arb_mul_ui(co,co,2,p);
    acb_mul_arb(v,v,co,p);
    acb_add(reg,reg,v,p);
    arb_mul(tx,tx,o->x,p);
    arb_mul(tx,tx,o->x,p);
  }
  acb_add(out,out,reg,p);
  for(int n=1;n<=32;n++){
    arb_mul_ui(nx,o->x,n,p);
    arb_sin_cos(sn,cn,nx,p);
    arb_mul(co,o->a,sn,p);
    arb_mul_ui(co,co,2,p);
    arb_div_ui(co,co,n,p);
    arb_sub(co,cn,co,p);
    arb_mul_2exp_si(nx,nx,-1);
    arb_sin(nx,nx,p);
    arb_sqr(nx,nx,p);
    arb_mul(nx,nx,o->d,p);
    arb_mul_ui(nx,nx,2,p);
    arb_div_ui(nx,nx,n*n,p);
    arb_add(co,co,nx,p);
    arb_mul(co,co,s[n],p);
    arb_mul_ui(co,co,2,p);
    acb_mul_ui(v,y,n,p);
    acb_cos(v,v,p);
    acb_mul_arb(v,v,co,p);
    acb_add(out,out,v,p);
  }
  /* Global regular-series remainder, propagated through original K/I/J
    combination; tiny though deliberately not sharp. */  acb_sub_arb(v,y,pi,p);
  acb_abs(rad,v,p);
  arb_add(rad,rad,o->x,p);
  tail(e,rad);
  tail(aa,pi);
  arb_add(e,e,aa,p);
  arb_add(aa,rad,pi,p);
  arb_abs(bb,o->a);
  arb_mul(bb,bb,aa,p);
  arb_mul_ui(bb,bb,4,p);
  arb_add_ui(bb,bb,2,p);
  arb_sqr(aa,aa,p);
  arb_abs(fac,o->d);
  arb_mul(fac,fac,aa,p);
  arb_mul_ui(fac,fac,8,p);
  arb_add(bb,bb,fac,p);
  arb_mul(e,e,bb,p);
  /* Smooth Fourier tail: bound brackets by 1+2|a|/n+2|d|/n². */  arb_abs(aa,acb_imagref(y));
  arb_set_ui(bb,2);
  arb_sub(aa,bb,aa,p);
  arb_neg(aa,aa);
  arb_exp(aa,aa,p);
  arb_one(bb);
  arb_sub(bb,bb,aa,p);
  if(!arb_is_positive(bb)){
    arb_indeterminate(aa);
  }
  else {
    arb_pow_ui(aa,aa,33,p);
    arb_div(aa,aa,bb,p);
  }
  arb_mul_ui(aa,aa,4,p);
  arb_div(aa,aa,pi,p);
  arb_abs(bb,o->a);
  arb_mul_ui(bb,bb,2,p);
  arb_abs(fac,o->d);
  arb_mul_ui(fac,fac,2,p);
  arb_add(bb,bb,fac,p);
  arb_add_ui(bb,bb,1,p);
  arb_mul(aa,aa,bb,p);
  arb_add(e,e,aa,p);
  acb_add_error_arb(out,e);
  CC(z);
  CC(r);
  CC(r2);
  CC(pow);
  CC(weight);
  CC(tmp);
  CC(v);
  CC(yy);
  CC(reg);
  AC(A);
  AC(B);
  AC(b);
  AC(co);
  AC(tx);
  AC(nx);
  AC(sn);
  AC(cn);
  AC(e);
  AC(rad);
  AC(aa);
  AC(bb);
  AC(rr);
  AC(fac);
}
/* Inner row with its homogeneous singular terms scaled before enclosure. */ static void inner_scaled(acb_t out,const acb_t t,CT*o){
  acb_t t2,w,ss,rs,ws,v,b0,b1,b2,tmp,z,reg;
  arb_t A,B,H,fac,xpow,co,e,rad,nx,sn,cn;
  CI(t2);
  CI(w);
  CI(ss);
  CI(rs);
  CI(ws);
  CI(v);
  CI(b0);
  CI(b1);
  CI(b2);
  CI(tmp);
  CI(z);
  CI(reg);
  AI(A);
  AI(B);
  AI(H);
  AI(fac);
  AI(xpow);
  AI(co);
  AI(e);
  AI(rad);
  AI(nx);
  AI(sn);
  AI(cn);
  arb_set(A,o->A);
  arb_set(B,o->B);
  arb_mul(H,o->H,c,p);
  acb_sqr(t2,t,p);
  acb_one(w);
  acb_sub(w,w,t2,p);
  acb_neg(ss,t2);
  acb_add_ui(ss,ss,2,p);
  acb_sqrt_analytic(rs,ss,analytic,p);
  acb_sqrt_analytic(ws,w,analytic,p);
  acb_inv(v,rs,p);
  acb_sub_ui(v,v,2,p);
  acb_add(b1,t,rs,p);
  acb_sub_ui(b1,b1,2,p);
  acb_mul_arb(b1,b1,A,p);
  acb_mul_ui(b1,b1,4,p);
  acb_sub(v,v,b1,p);
  acb_mul(b2,t,t2,p);
  acb_mul(tmp,ss,rs,p);
  acb_add(b2,b2,tmp,p);
  acb_sub_ui(b2,b2,2,p);
  acb_mul(tmp,w,ws,p);
  acb_mul_ui(tmp,tmp,2,p);
  acb_sub(b2,b2,tmp,p);
  acb_mul_arb(b2,b2,B,p);
  acb_mul_ui(b2,b2,4,p);
  acb_div_ui(b2,b2,3,p);
  acb_add(v,v,b2,p);
  acb_mul(v,v,t,p);
  acb_add_ui(v,v,1,p);
  acb_mul_ui(v,v,2,p);
  acb_div(v,v,w,p);
  acb_mul_arb(out,v,H,p);
  if(o->op==1){
    acb_sub(v,t,rs,p);
    acb_mul_ui(tmp,ws,2,p);
    acb_add(v,v,tmp,p);
    acb_mul_arb(v,v,A,p);
    acb_mul_ui(v,v,2,p);
    acb_mul(v,v,t,p);
    acb_div(tmp,t,rs,p);
    acb_add(v,v,tmp,p);
    acb_neg(v,v);
    acb_add_ui(v,v,1,p);
    acb_mul_ui(v,v,2,p);
    acb_mul_arb(out,v,H,p);
  }
  /* Inner regular power series. Differences depend only on t, not x. */  arb_one(xpow);
  acb_zero(reg);
  for(int m=1;m<M;m++){
    arb_mul(xpow,xpow,o->x,p);
    arb_mul(xpow,xpow,o->x,p);
    acb_pow_ui(b0,t2,2*m,p);
    acb_pow_ui(tmp,ss,2*m,p);
    acb_add(b0,b0,tmp,p);
    acb_sub_ui(b0,b0,2,p);
    acb_pow_ui(b1,t2,2*m+1,p);
    acb_pow_ui(tmp,ss,2*m+1,p);
    acb_add(b1,b1,tmp,p);
    acb_sub_ui(b1,b1,2,p);
    acb_mul_arb(b1,b1,A,p);
    acb_mul_ui(b1,b1,2,p);
    acb_div_ui(b1,b1,2*m+1,p);
    acb_sub(b0,b0,b1,p);
    acb_pow_ui(b2,t2,2*m+2,p);
    acb_pow_ui(tmp,ss,2*m+2,p);
    acb_add(b2,b2,tmp,p);
    acb_sub_ui(b2,b2,2,p);
    acb_pow_ui(tmp,w,2*m+2,p);
    acb_mul_ui(tmp,tmp,2,p);
    acb_sub(b2,b2,tmp,p);
    acb_mul_arb(b2,b2,B,p);
    acb_div_ui(b2,b2,(2*m+1)*(2*m+2),p);
    acb_add(b0,b0,b2,p);
    if(o->op==1){
      acb_pow_ui(b0,t2,2*m,p);
      acb_pow_ui(tmp,ss,2*m,p);
      acb_sub(b0,b0,tmp,p);
      acb_pow_ui(b1,t2,2*m+1,p);
      acb_pow_ui(tmp,ss,2*m+1,p);
      acb_sub(b1,b1,tmp,p);
      acb_pow_ui(tmp,w,2*m+1,p);
      acb_mul_ui(tmp,tmp,2,p);
      acb_add(b1,b1,tmp,p);
      acb_mul_arb(b1,b1,A,p);
      acb_div_ui(b1,b1,2*m+1,p);
      acb_sub(b0,b0,b1,p);
    }
    arb_mul(co,xpow,kl[m],p);
    acb_mul_arb(b0,b0,co,p);
    acb_add(reg,reg,b0,p);
  }
  for(int n=1;n<=32;n++){
    arb_mul_ui(nx,o->x,n,p);
    arb_sin_cos(sn,cn,nx,p);
    arb_mul(co,o->a,sn,p);
    arb_mul_ui(co,co,2,p);
    arb_div_ui(co,co,n,p);
    arb_sub(co,cn,co,p);
    arb_mul_2exp_si(nx,nx,-1);
    arb_sin(nx,nx,p);
    arb_sqr(nx,nx,p);
    arb_mul(nx,nx,o->d,p);
    arb_mul_ui(nx,nx,2,p);
    arb_div_ui(nx,nx,n*n,p);
    arb_add(co,co,nx,p);
    if(o->op==1){
      arb_mul_ui(nx,o->x,n,p);
      arb_sin(sn,nx,p);
      arb_mul_2exp_si(nx,nx,-1);
      arb_sin(nx,nx,p);
      arb_sqr(nx,nx,p);
      arb_mul(nx,nx,o->a,p);
      arb_mul_ui(nx,nx,2,p);
      arb_div_ui(nx,nx,n,p);
      arb_sub(co,sn,nx,p);
      arb_mul(co,co,s[n],p);
      arb_mul_ui(co,co,2,p);
      acb_mul_arb(v,w,o->x,p);
      acb_mul_ui(v,v,n,p);
      acb_sin(v,v,p);
    }
    else {
      arb_mul(co,co,s[n],p);
      arb_mul_si(co,co,-4,p);
      acb_mul_arb(v,w,o->x,p);
      acb_mul_ui(v,v,n,p);
      acb_mul_2exp_si(v,v,-1);
      acb_sin(v,v,p);
      acb_sqr(v,v,p);
    }
    acb_mul_arb(v,v,co,p);
    acb_add(reg,reg,v,p);
  }
  /* Clausen tail for |z| up to max(|x(1+-w)|,|xw|,x), using the
    original K/I/J differences as a conservative propagation envelope. */  acb_abs(rad,ss,p);
  acb_abs(e,t2,p);
  arb_max(rad,rad,e,p);
  acb_abs(e,w,p);
  arb_max(rad,rad,e,p);
  arb_one(e);
  arb_max(rad,rad,e,p);
  arb_mul(rad,rad,o->x,p);
  arb_div(e,rad,L,p);
  arb_sqr(e,e,p);
  arb_one(fac);
  arb_sub(fac,fac,e,p);
  if(!arb_is_positive(fac))arb_indeterminate(e);
  else{
    arb_pow_ui(e,e,M,p);
    arb_mul(e,e,c,p);
    arb_mul_ui(e,e,4,p);
    arb_div(e,e,rootL,p);
    arb_div(e,e,fac,p);
  }
  arb_abs(fac,o->a);
  arb_mul(fac,fac,rad,p);
  arb_mul_ui(fac,fac,8,p);
  arb_add_ui(fac,fac,4,p);
  arb_sqr(co,rad,p);
  arb_abs(nx,o->d);
  arb_mul(co,co,nx,p);
  arb_mul_ui(co,co,6,p);
  arb_add(fac,fac,co,p);
  arb_mul(e,e,fac,p);
  /* complex Fourier tail on xw and x(1+-w). */  acb_mul_arb(v,w,o->x,p);
  arb_abs(co,acb_imagref(v));
  arb_set_ui(nx,2);
  arb_sub(co,nx,co,p);
  if(!arb_is_positive(co))arb_indeterminate(e);
  else{
    arb_neg(co,co);
    arb_exp(co,co,p);
    arb_one(nx);
    arb_sub(nx,nx,co,p);
    arb_pow_ui(co,co,33,p);
    arb_div(co,co,nx,p);
    arb_mul_ui(co,co,8,p);
    arb_div(co,co,pi,p);
    arb_abs(nx,o->a);
    arb_mul_ui(nx,nx,2,p);
    arb_abs(sn,o->d);
    arb_mul_ui(sn,sn,2,p);
    arb_add(nx,nx,sn,p);
    arb_add_ui(nx,nx,1,p);
    arb_mul(co,co,nx,p);
    arb_add(e,e,co,p);
  }
  acb_add_error_arb(reg,e);
  acb_mul(reg,reg,t,p);
  if(o->op==2)acb_div(reg,reg,w,p);
  acb_mul_arb(reg,reg,o->x,p);
  acb_div_arb(reg,reg,o->u,p);
  acb_add(out,out,reg,p);
  CC(t2);
  CC(w);
  CC(ss);
  CC(rs);
  CC(ws);
  CC(v);
  CC(b0);
  CC(b1);
  CC(b2);
  CC(tmp);
  CC(z);
  CC(reg);
  AC(A);
  AC(B);
  AC(H);
  AC(fac);
  AC(xpow);
  AC(co);
  AC(e);
  AC(rad);
  AC(nx);
  AC(sn);
  AC(cn);
}
static int fun(acb_ptr out,const acb_t t,void*param,slong order,slong prec){
  (void)prec;
  analytic=order;
  CT*o=param;
  calls++;
  acb_t ts,y,zm,zp,r0,rL,k1,i1,j1,k2,i2,j2,ky,iy,jy,v,jac,z;
  arb_t root;
  CI(ts);
  CI(y);
  CI(zm);
  CI(zp);
  CI(r0);
  CI(rL);
  CI(k1);
  CI(i1);
  CI(j1);
  CI(k2);
  CI(i2);
  CI(j2);
  CI(ky);
  CI(iy);
  CI(jy);
  CI(v);
  CI(jac);
  CI(z);
  AI(root);
  if(!o->region){
    arb_set_str(root,"1.5",p);
    if(arb_lt(o->x,root)){
      inner_scaled(out,t,o);
      goto signfinish;
    }
  }
  arb_sqrt(root,o->w,p);
  acb_sqr(ts,t,p);
  acb_mul_arb(zm,ts,o->w,p);
  acb_mul_arb(r0,t,root,p);
  acb_mul_arb(jac,t,o->w,p);
  acb_mul_ui(jac,jac,2,p);
  if(!o->region){
    acb_one(y);
    acb_sub(y,y,ts,p);
    acb_mul_arb(y,y,o->x,p);
  }
  else acb_add_arb(y,zm,o->x,p);
  acb_add_arb(zp,y,o->x,p);
  if(o->region&&o->op==2){
    arb_t ay,rx;
    AI(ay);
    AI(rx);
    acb_abs(ay,y,p);
    arb_div(rx,o->x,ay,p);
    arb_set_str(root,".65",p);
    int stable=arb_lt(rx,root);
    AC(ay);
    AC(rx);
    arb_sqrt(root,o->w,p);
    if(stable){
      outer_raw(out,y,o);
      acb_mul(out,out,jac,p);
      acb_div(out,out,y,p);
      acb_mul_arb(out,out,o->x,p);
      acb_div_arb(out,out,o->u,p);
      acb_div_ui(out,out,2,p);
      goto signfinish;
    }
  }
  kernel(k1,i1,j1,zm,r0,NULL,1);
  if(o->region){
    acb_neg(i1,i1);
    acb_neg(z,ts);
    acb_add_ui(z,z,2,p);
    acb_sqrt_analytic(rL,z,analytic,p);
    acb_mul_arb(rL,rL,root,p);
    kernel(k2,i2,j2,zp,NULL,rL,2);
  }
  else kernel(k2,i2,j2,zp,NULL,NULL,2);
  kernel(ky,iy,jy,y,NULL,NULL,0);
  if(o->op==1){
    acb_sub(out,k1,k2,p);
    acb_sub(v,i1,i2,p);
    acb_mul_ui(z,iy,2,p);
    acb_add(v,v,z,p);
    acb_mul_arb(v,v,o->a,p);
    acb_sub(out,out,v,p);
  }
  else{
    acb_add(out,k1,k2,p);
    acb_add(v,i1,i2,p);
    acb_mul_arb(v,v,o->a,p);
    acb_mul_ui(v,v,2,p);
    acb_sub(out,out,v,p);
    acb_add(v,j1,j2,p);
    acb_mul_ui(z,jy,2,p);
    acb_sub(v,v,z,p);
    acb_mul_arb(v,v,o->d,p);
    acb_add(out,out,v,p);
    if(!o->region)acb_sub_arb(out,out,o->f,p);
  }
  acb_mul(out,out,jac,p);
  arb_mul(root,root,c,p);
  arb_mul_ui(root,root,2,p);
  acb_add_arb(out,out,root,p);
  if(o->region){
    acb_neg(z,ts);
    acb_add_ui(z,z,2,p);
    acb_sqrt_analytic(z,z,analytic,p);
    acb_div(v,t,z,p);
    acb_mul_arb(v,v,root,p);
    if(o->op==1)acb_sub(out,out,v,p);
    else acb_add(out,out,v,p);
  }
  else{
    acb_neg(z,zp);
    acb_add_arb(z,z,L,p);
    if(!order&&arb_contains_zero(acb_realref(z))&&arb_is_zero(acb_imagref(t))){
      /* 0<=t/sqrt(2pi-2x+x t²)<=1/sqrt(x), uniformly up to x=pi. */     arb_t lo,hi;
      AI(lo);
      AI(hi);
      arb_zero(lo);
      arb_sqrt(hi,o->x,p);
      arb_inv(hi,hi,p);
      if(!arb_is_nonnegative(acb_realref(t)))arb_neg(lo,hi);
      arb_union(acb_realref(v),lo,hi,p);
      arb_zero(acb_imagref(v));
      AC(lo);
      AC(hi);
    }
    else{
      acb_sqrt_analytic(z,z,analytic,p);
      acb_div(v,t,z,p);
    }
    acb_mul_arb(v,v,o->x,p);
    acb_mul_arb(v,v,c,p);
    acb_mul_ui(v,v,2,p);
    if(o->op==1)acb_sub(out,out,v,p);
    else acb_add(out,out,v,p);
  }
  if(o->op==2){
    acb_div(out,out,y,p);
    acb_mul_arb(out,out,o->x,p);
  }
  acb_div_arb(out,out,o->u,p);
  acb_div_ui(out,out,2,p);
  signfinish:  if(order){
    if(arb_is_negative(acb_realref(out)))acb_neg(out,out);
    else if(!arb_is_positive(acb_realref(out)))acb_indeterminate(out);
  }
  else{
    arb_abs(acb_realref(out),acb_realref(out));
    arb_zero(acb_imagref(out));
  }
  CC(ts);
  CC(y);
  CC(zm);
  CC(zp);
  CC(r0);
  CC(rL);
  CC(k1);
  CC(i1);
  CC(j1);
  CC(k2);
  CC(i2);
  CC(j2);
  CC(ky);
  CC(iy);
  CC(jy);
  CC(v);
  CC(jac);
  CC(z);
  AC(root);
  return 0;
}
static void initct(CT*o,const arb_t x,int op){
  AI(o->x);
  AI(o->u);
  AI(o->a);
  AI(o->d);
  AI(o->f);
  AI(o->w);
  AI(o->A);
  AI(o->B);
  AI(o->H);
  o->op=op;
  arb_set(o->x,x);
  arb_t du,dd,t;
  AI(du);
  AI(dd);
  AI(t);
  compute_u0(o->u,x,&up,p);
  compute_du0(du,x,&up,p);
  compute_ddu0(dd,x,&up,p);
  arb_div(o->a,du,o->u,p);
  arb_sqr(o->d,o->a,p);
  arb_mul_ui(o->d,o->d,2,p);
  arb_div(t,dd,o->u,p);
  arb_sub(o->d,o->d,t,p);
  Kp(o->f,x,p);
  IKp(t,x,p);
  arb_mul(t,t,o->a,p);
  arb_mul_ui(t,t,2,p);
  arb_sub(o->f,o->f,t,p);
  IIKp(t,x,p);
  arb_mul(t,t,o->d,p);
  arb_add(o->f,o->f,t,p);
  arb_mul_ui(o->f,o->f,2,p);
  arb_mul(o->A,x,o->a,p);
  arb_sqr(o->B,x,p);
  arb_mul(o->B,o->B,o->d,p);
  arb_sqrt(o->H,x,p);
  arb_div(o->H,o->H,o->u,p);
  arb_mul_2exp_si(o->H,o->H,-1);
  arb_ptr jets=_arb_vec_init(3);
  if(derivative_profile_jets(jets,x,&up,p)){
    arb_div(o->A,jets+1,jets,p);
    arb_sqr(o->B,o->A,p);
    arb_mul_ui(o->B,o->B,2,p);
    arb_div(t,jets+2,jets,p);
    arb_sub(o->B,o->B,t,p);
    arb_inv(o->H,jets,p);
    arb_mul_2exp_si(o->H,o->H,-1);
  }
  _arb_vec_clear(jets,3);
  AC(du);
  AC(dd);
  AC(t);
}
static void clearct(CT*o){
  AC(o->x);
  AC(o->u);
  AC(o->a);
  AC(o->d);
  AC(o->f);
  AC(o->w);
  AC(o->A);
  AC(o->B);
  AC(o->H);
}
static void nearzero(arb_t out,CT*o,const arb_t eta){
  arb_t h,z,t,b,v,r,e,lo,hi;
  AI(h);
  AI(z);
  AI(t);
  AI(b);
  AI(v);
  AI(r);
  AI(e);
  AI(lo);
  AI(hi);
  arb_mul(h,o->x,eta,p);
  arb_sub(lo,o->x,h,p);
  arb_add(hi,o->x,h,p);
  arb_union(z,lo,hi,p);
  ddKp(b,z,p);
  arb_abs(b,b);
  dKp(v,z,p);
  arb_abs(v,v);
  arb_abs(t,o->a);
  arb_mul(v,v,t,p);
  arb_mul_ui(v,v,2,p);
  arb_add(b,b,v,p);
  Kp(v,z,p);
  arb_abs(v,v);
  arb_abs(t,o->d);
  arb_mul(v,v,t,p);
  arb_add(b,b,v,p);
  arb_sqr(z,h,p);
  arb_mul(b,b,z,p);
  arb_div_ui(b,b,2,p);
  arb_zero(lo);
  arb_union(hi,lo,h,p);
  acb_t az,K,I,J;
  CI(az);
  CI(K);
  CI(I);
  CI(J);
  acb_set_arb(az,hi);
  kernel(K,I,J,az,NULL,NULL,1);
  arb_abs(r,acb_realref(K));
  arb_mul(v,r,z,p);
  arb_div_ui(v,v,4,p);
  arb_sqrt(e,h,p);
  arb_mul(e,e,h,p);
  arb_mul(e,e,c,p);
  arb_mul_ui(e,e,8,p);
  arb_div_ui(e,e,9,p);
  arb_add(v,v,e,p);
  arb_mul(v,v,t,p);
  arb_mul_ui(v,v,2,p);
  arb_add(b,b,v,p);
  arb_mul(b,b,o->x,p);
  arb_div(b,b,o->u,p);
  arb_div_ui(out,b,2,p);
  CC(az);
  CC(K);
  CC(I);
  CC(J);
  AC(h);
  AC(z);
  AC(t);
  AC(b);
  AC(v);
  AC(r);
  AC(e);
  AC(lo);
  AC(hi);
}
static void diagonal_tail(arb_t out,CT*o,const arb_t tau){
  arb_t h,t,b,v;
  AI(h);
  AI(t);
  AI(b);
  AI(v);
  arb_sqr(h,tau,p);
  arb_mul(h,h,o->x,p);
  IKp(out,h,p);
  arb_mul_ui(out,out,2,p);
  arb_abs(b,o->a);
  arb_mul_ui(b,b,4,p);
  if(o->op==2){
    arb_abs(t,o->d);
    arb_mul(t,t,pi,p);
    arb_mul_ui(t,t,4,p);
    arb_add(b,b,t,p);
    arb_abs(t,o->f);
    arb_add(b,b,t,p);
  }
  arb_mul(b,b,h,p);
  arb_add(out,out,b,p);
  arb_div(out,out,o->u,p);
  arb_div_ui(out,out,2,p);
  if(o->op==2){
    arb_sqr(t,tau,p);
    arb_one(v);
    arb_sub(v,v,t,p);
    arb_div(out,out,v,p);
  }
  AC(h);
  AC(t);
  AC(b);
  AC(v);
}
static int row(arb_t out,const arb_t x,int op,int verbose){
  CT o;
  initct(&o,x,op);
  acb_t aa,bb,res;
  CI(aa);
  CI(bb);
  CI(res);
  arb_t eta,e,t;
  AI(eta);
  AI(e);
  AI(t);
  arb_set_str(eta,".005",p);
  mag_t tol;
  mag_init(tol);
  int goal=op==1?16:12;
  double tol_d=op==1?1e-5:1e-4;
  if(arf_cmp_d(arb_midref(x),1.4)>0){
    goal=9;
    tol_d=.0005;
  }
  if(op==1&&arf_cmp_d(arb_midref(x),3.1)>0){
    goal=6;
    tol_d=.01;
  }
  mag_set_d(tol,tol_d);
  acb_calc_integrate_opt_t opt;
  acb_calc_integrate_opt_init(opt);
  opt->eval_limit=3000;
  opt->depth_limit=25;
  opt->deg_limit=40;
  arb_zero(out);
  int ok=arb_is_positive(o.u);
  for(int r=0;r<2&&ok;r++){
    o.region=r;
    if(r&&op==1){
      /* Outer positivity follows modewise from x*u'/u<=2. Use the
      inherited exact primitive formula, audited independently. */    arb_set(t,o.A);
      arb_set_ui(e,2);
      if(!arb_le(t,e)){
        ok=0;
        break;
      }
      arb_set_str(t,"3.1",p);
      if(arb_gt(x,t)){
        arb_sub(t,pi,x,p);
        arb_nonnegative_part(t,t);
        arb_get_ubound_arf(arb_midref(t),t,p);
        mag_zero(arb_radref(t));
        IKp(e,t,p);
        arb_mul_ui(e,e,2,p);
        arb_abs(o.w,o.a);
        arb_mul(o.w,o.w,t,p);
        arb_mul_ui(o.w,o.w,4,p);
        arb_add(e,e,o.w,p);
        arb_div(e,e,o.u,p);
        arb_div_ui(e,e,2,p);
        arb_add(out,out,e,p);
        ok=arb_is_finite(out);
        continue;
      }
      Integration_params_struct ip;
      ip.up=up;
      arb_init(ip.x);
      arb_set(ip.x,x);
      ip.order=2;
      ip.u=_arb_vec_init(3);
      arb_set(ip.u,o.u);
      arb_mul(ip.u+1,o.a,o.u,p);
      compute_ddu0(ip.u+2,x,&up,p);
      compute_big_y_T1(e,&ip,p);
      arb_add(out,out,e,p);
      if(verbose){
        printf("OP1 exact outer ");
        arb_printn(e,16,0);
        puts("");
      }
      ok=arb_is_finite(out);
      arb_clear(ip.x);
      _arb_vec_clear(ip.u,3);
      continue;
    }
    if(r)arb_sub(o.w,pi,x,p);
    else arb_set(o.w,x);
    if(arb_is_zero(o.w))continue;
    if(r)arb_nonnegative_part(o.w,o.w);
    acb_zero(aa);
    acb_one(bb);
    arb_set_str(t,"1.4",p);
    if(!r&&arb_gt(x,t)){
      arb_set_str(t,".005",p);
      acb_set_arb(aa,t);
      diagonal_tail(e,&o,t);
      arb_add(out,out,e,p);
    }
    if(!r&&op==2){
      arb_one(t);
      arb_sub(t,t,eta,p);
      arb_sqrt(t,t,p);
      acb_set_arb(bb,t);
      nearzero(e,&o,eta);
      arb_add(out,out,e,p);
    }
    int status=acb_calc_integrate(res,fun,&o,aa,bb,goal,tol,opt,p);
    if(verbose){
      printf("OP%d region%d status%d calls%ld ",op,r,status,calls);
      acb_printn(res,14,0);
      puts("");
      fflush(stdout);
    }
    arb_add(out,out,acb_realref(res),p);
    ok=arb_is_finite(out)&&acb_is_finite(res)&&status==0;
  }
  clearct(&o);
  CC(aa);
  CC(bb);
  CC(res);
  AC(eta);
  AC(e);
  AC(t);
  mag_clear(tol);
  return ok;
}
static int cover(const arb_t l,const arb_t r,int op,int depth){
  arb_t x,m,v;
  AI(x);
  AI(m);
  AI(v);
  arb_union(x,l,r,p);
  nodes++;
  arb_mul_ui(m,l,103,p);
  arb_div_ui(m,m,100,p);
  int narrow=arb_le(r,m);
  int ok=narrow&&row(v,x,op,0)&&arb_lt(v,lim);
  if(ok){
    leaves++;
    printf("OP%d LEAF %ld x=",op,leaves);
    arb_printn(x,12,0);
    printf(" row=");
    arb_printn(v,12,0);
    printf(" calls=%ld\n",calls);
    fflush(stdout);
  }
  else if(depth<30&&nodes<100000){
    arb_add(m,l,r,p);
    arb_mul_2exp_si(m,m,-1);
    ok=cover(l,m,op,depth+1);
    if(ok)ok=cover(m,r,op,depth+1);
  }
  else{
    printf("UNRESOLVED OP%d ",op);
    arb_printn(x,16,0);
    printf(" ");
    arb_printn(v,16,0);
    puts("");
  }
  AC(x);
  AC(m);
  AC(v);
  return ok;
}
int main(int argc,char**argv){
  setup();
  arb_t x,y;
  AI(x);
  AI(y);
  if(argc>1&&!strcmp(argv[1],"crosscheck")){
    int ok=1;
    for(int j=1;j<=60;j++){
      arb_mul_ui(x,pi,j,p);
      arb_div_ui(x,x,31,p);
      acb_t z,K,I,J;
      CI(z);
      CI(K);
      CI(I);
      CI(J);
      acb_set_arb(z,x);
      kernel(K,I,J,z,NULL,NULL,0);
      Kp(y,x,p);
      ok=ok&&arb_overlaps(y,acb_realref(K));
      IKp(y,x,p);
      ok=ok&&arb_overlaps(y,acb_realref(I));
      IIKp(y,x,p);
      ok=ok&&arb_overlaps(y,acb_realref(J));
      CC(z);
      CC(K);
      CC(I);
      CC(J);
    }
    printf("global kernel crosscheck %s\n",ok?"PASS":"FAIL");
    return ok?0:2;
  }
  int op=argc>1?atoi(argv[1]):2;
  if(argc>2&&!strcmp(argv[2],"pi"))arb_set(x,pi);
  else arb_set_str(x,argc>2?argv[2]:".1",p);
  if(argc>3&&!strcmp(argv[3],"reject")){
    int ok=row(y,x,op,0);
    arb_set_str(lim,".1",p);
    int accepted=ok&&arb_lt(y,lim);
    printf("IMPOSSIBLE_TARGET %s\n",accepted?"INVALID_ACCEPT":"UNRESOLVED");
    return accepted?0:2;
  }
  if(argc>3&&!strcmp(argv[3],"cover")){
    arb_set_str(lim,argc>4?argv[4]:".9974",p);
    int ok=cover(x,pi,op,0);
    printf("OP%d COMPACT %s leaves=%ld nodes=%ld calls=%ld\n",op,ok?"PASS":"UNRESOLVED",leaves,nodes,calls);
    return ok?0:2;
  }
  if(argc>3&&!strcmp(argv[3],"sample")){
    CT o;
    initct(&o,x,op);
    o.region=0;
    arb_set(o.w,x);
    acb_t tt,res;
    CI(tt);
    CI(res);
    for(int j=0;j<11;j++){
      acb_set_ui(tt,j);
      acb_div_ui(tt,tt,10,p);
      fun(res,tt,&o,0,p);
      acb_printn(tt,10,0);
      printf(" ");
      acb_printn(res,20,0);
      puts("");
    }
    return 0;
  }
  if(argc>3){
    arb_set_str(y,argv[3],p);
    arb_add_error(x,y);
  }
  int ok=row(y,x,op,1);
  printf("OP%d POINT %s ",op,ok?"UPPER_BOUND":"UNRESOLVED");
  arb_printn(y,20,0);
  puts("");
  arb_get_ubound_arf(arb_midref(x),y,p);
  mag_zero(arb_radref(x));
  printf("ROW_UPPER ");
  arb_printn(x,25,0);
  puts("");
  return ok?0:2;
}
