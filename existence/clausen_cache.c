/* AI-added acceleration of the inherited Clausen evaluator. See PROOF.md.
 * Fixed-precision, single-thread cache; no values from samples enter a proof. */
#include "basic_functions.h"
#define CM 48
#define NC 128
static struct {
  int used, sine;
  slong prec;
  arb_t s,g,tail;
  arb_ptr a;
}
cache[NC];
int clausen_cached(arb_t out,const arb_t s,const arb_t input,int sine,slong p) {
  if (arb_is_int(s)) return 0;
  /* inherited direct evaluator handles integer limits */
  arb_t x,pi,two,t,z,r,d,e,f,sg;
  arb_init(x);
  arb_init(pi);
  arb_init(two);
  arb_init(t);
  arb_init(z);
  arb_init(r);
  arb_init(d);
  arb_init(e);
  arb_init(f);
  arb_init(sg);
  arb_const_pi(pi,p);
  arb_mul_2exp_si(two,pi,1);
  arb_abs(x,input);
  arb_one(sg);
  if (sine && arb_is_negative(input)) arb_neg(sg,sg);
  arb_set_si(t,-2);
  arb_set_ui(z,7);
  int ok=arb_ge(s,t)&&arb_le(s,z)&&arb_le(x,z);
  if(!ok) goto done;
  /* Periodic reflection. Splitting is unnecessary: absolute value encloses
    the folded coordinate; sine signs are enclosed when crossing pi. */
  if(arb_gt(x,pi)){
    arb_sub(x,two,x,p);
    if(sine&&!arb_is_nonnegative(x))arb_zero_pm_one(sg);
    arb_abs(x,x);
    arb_neg(sg,sg);
  }
  else if(!arb_lt(x,pi)){
    arb_sub(x,x,pi,p);
    arb_abs(x,x);
    arb_sub(x,pi,x,p);
    if(sine)arb_zero_pm_one(sg);
  }
  arb_nonnegative_part(x,x);
  /* The rounded enclosure of pi may extend above pi. The same
   * convergent tail formula is valid on [0,3.2], with R<1. */
  arb_set_str(t,"3.2",p);
  if(!arb_le(x,t)){
    ok=0;
    goto done;
  }
  if(!arb_is_positive(x)){
    arb_one(t);
    if(!arb_gt(s,t)){
      ok=0;
      goto done;
    }
  }
  int n;
  for(n=0;n<NC;n++) if(!cache[n].used || (cache[n].prec==p&&cache[n].sine==sine&&arb_equal(cache[n].s,s)))break;
  if(n==NC){
    ok=0;
    goto done;
  }
  if(!cache[n].used){
    cache[n].used=1;
    cache[n].prec=p;
    cache[n].sine=sine;
    arb_init(cache[n].s);
    arb_init(cache[n].g);
    arb_init(cache[n].tail);
    cache[n].a=_arb_vec_init(CM);
    arb_set(cache[n].s,s);
    arb_neg(t,s);
    arb_add_ui(t,t,1,p);
    arb_gamma(cache[n].g,t,p);
    arb_mul_2exp_si(t,s,-1);
    if(sine)arb_cos_pi(t,t,p);
    else arb_sin_pi(t,t,p);
    arb_mul(cache[n].g,cache[n].g,t,p);
    for(int m=0;m<CM;m++){
      arb_sub_ui(t,s,2*m+sine,p);
      arb_zeta(cache[n].a+m,t,p);
      arb_fac_ui(t,2*m+sine,p);
      arb_div(cache[n].a+m,cache[n].a+m,t,p);
      if(m%2)arb_neg(cache[n].a+m,cache[n].a+m);
    }
    /* zeta decreases above 1; Gamma ratio <= (k+2)^2 for -2<=s<=7,
     k>=96. The functional equation gives the following common prefactor. */
    arb_neg(t,s);
    arb_add_ui(t,t,2*CM+1+sine,p);
    arb_zeta(cache[n].tail,t,p);
    arb_sub_ui(t,s,1,p);
    arb_pow(t,two,t,p);
    arb_mul(cache[n].tail,cache[n].tail,t,p);
    arb_mul_ui(cache[n].tail,cache[n].tail,2,p);
  }
  arb_sqr(z,x,p);
  arb_set(out,cache[n].a+CM-1);
  for(int m=CM-2;m>=0;m--){
    arb_mul(out,out,z,p);
    arb_add(out,out,cache[n].a+m,p);
  }
  if(sine)arb_mul(out,out,x,p);
  arb_sub_ui(t,s,1,p);
  arb_pow(t,x,t,p);
  arb_mul(t,t,cache[n].g,p);
  arb_add(out,out,t,p);
  arb_div(r,x,two,p);
  arb_sqr(r,r,p);
  arb_neg(d,r);
  arb_add_ui(d,d,1,p);
  /* sum_{m>=M}(2m+2+sine)^2 r^m, multiplied by (x/(2pi))^sine */
  int b=2*CM+2+sine;
  arb_set_ui(e,b*b);
  arb_div(e,e,d,p);
  arb_sqr(t,d,p);
  arb_mul_ui(f,r,4*b,p);
  arb_div(f,f,t,p);
  arb_add(e,e,f,p);
  arb_mul(t,t,d,p);
  arb_add_ui(f,r,1,p);
  arb_mul(f,f,r,p);
  arb_mul_ui(f,f,4,p);
  arb_div(f,f,t,p);
  arb_add(e,e,f,p);
  arb_pow_ui(t,r,CM,p);
  arb_mul(e,e,t,p);
  arb_mul(e,e,cache[n].tail,p);
  if(sine){
    arb_div(t,x,two,p);
    arb_mul(e,e,t,p);
  }
  arb_add_error(out,e);
  if(sine)arb_mul(out,out,sg,p);
  ok=arb_is_finite(out);
  done:arb_clear(x);
  arb_clear(pi);
  arb_clear(two);
  arb_clear(t);
  arb_clear(z);
  arb_clear(r);
  arb_clear(d);
  arb_clear(e);
  arb_clear(f);
  arb_clear(sg);
  return ok;
}
