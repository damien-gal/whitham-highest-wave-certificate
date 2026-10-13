/* complex_functions.c: pointwise evaluators built on u_0 and the kernel.
 * See complex_functions.h for the list of routines.
 */
#include "complex_functions.h"


int compute_alpha1(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0,du0;
  arb_init(u0); arb_init(du0);
  
  compute_u0(u0,inp,params,prec);   
  compute_du0(du0,inp,params,prec);
    
  arb_sqr(out,inp,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,u0,prec);
  arb_div(out,out,u0,prec);
  arb_mul(out,out,du0,prec);
  arb_abs(out,out);
  
  arb_clear(u0);
  arb_clear(du0);
  
  return 0;
}

int compute_alpha2(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0,du0,ddu0;
  arb_t factor1,factor2;
  arb_init(u0); arb_init(du0); arb_init(ddu0);
  arb_init(factor1); arb_init(factor2);
  
  compute_u0(u0,inp,params,prec);   
  compute_du0(du0,inp,params,prec); 
  compute_ddu0(ddu0,inp,params,prec);
  
  // factor1 = x^3/(2*u0^2)
  arb_div(factor1,inp,u0,prec);
  arb_sqr(factor1,factor1,prec);
  arb_mul(factor1,factor1,inp,prec);
  arb_div_si(factor1,factor1,2,prec);
  
  // factor2 = u0'' - 2*(u0')^2/u0
  arb_sqr(factor2,du0,prec);
  arb_mul_si(factor2,factor2,2,prec);
  arb_div(factor2,factor2,u0,prec);
  arb_sub(factor2,ddu0,factor2,prec);
  
  arb_mul(out,factor1,factor2,prec);
  arb_abs(out,out);
  
  arb_clear(u0);
  arb_clear(du0);
  arb_clear(ddu0);
  
  arb_clear(factor1);
  arb_clear(factor2);
  
  return 0;
}

// f(x) = 2K + 2u'/u*K2(x,0) + 1/u^2*(2*(u')^2-u*u'')*K2bar(x,0)
// K2(x,0) = -2*IK(x)
// K2bar(x,0) = 2*IIK(x)
// f(x) = 2(K - 2u'/u*IK(x) + 1/u^2*(2*(u')^2-u*u'')*IIK(x))
int compute_alphaf(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0, du0, ddu0;
  arb_t aux,aux2;
  arb_init(u0); arb_init(du0); arb_init(ddu0);
  arb_init(aux); arb_init(aux2);
  
  compute_u0(u0,inp,params,prec);
  compute_du0(du0,inp,params,prec);
  compute_ddu0(ddu0,inp,params,prec);
  
  arb_mul(aux,du0,du0,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,u0,ddu0,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,u0,prec);
  arb_div(aux,aux,u0,prec);
  IIKp(aux2,inp,prec);
  
  arb_mul(out,aux,aux2,prec);
 

  arb_div(aux,du0,u0,prec);
  IKp(aux2,inp,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  Kp(aux,inp,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(out,out,inp,prec);
  arb_div(out,out,u0,prec);
  
  arb_abs(out,out);
  
  arb_clear(u0); arb_clear(du0); arb_clear(ddu0);
  arb_clear(aux); arb_clear(aux2);
  return 0;
}


int compute_convexity(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t ddu0;
  arb_init(ddu0);
  
  compute_ddu0(ddu0,inp,params,prec);
  arb_mul(out,ddu0,inp,prec);
  
  arb_clear(ddu0);
  
  return 0;
}


int compute_F0(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0,Lu0;
  arb_init(u0); arb_init(Lu0);
  
  compute_u0(u0,inp,params,prec);   
  compute_Lu0(Lu0,inp,params,prec);
  
  arb_div(out,Lu0,u0,prec);
  arb_sub(out,out,u0,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,inp,prec);
    
  arb_clear(u0);
  arb_clear(Lu0);
  
  return 0;
}


int compute_dF0(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0,du0;
  arb_t Lu0,dLu0;
  arb_t aux;
  arb_init(u0); arb_init(du0);
  arb_init(Lu0); arb_init(dLu0);
  arb_init(aux);
  
  compute_u0(u0,inp,params,prec); compute_du0(du0,inp,params,prec);  
  compute_Lu0(Lu0,inp,params,prec); compute_dLu0(dLu0,inp,params,prec);
  
  arb_div(aux,u0,inp,prec);
  arb_sub(out,aux,du0,prec);
  
  arb_div(aux,dLu0,u0,prec);
  arb_add(out,out,aux,prec);
  
  arb_div(aux,Lu0,u0,prec);
  arb_div(aux,aux,inp,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(aux,Lu0,du0,prec);
  arb_div(aux,aux,u0,prec);
  arb_div(aux,aux,u0,prec);
  arb_sub(out,out,aux,prec);
  
  arb_div(out,out,inp,prec);  
  arb_div_si(out,out,2,prec);
    
  arb_clear(u0); arb_clear(du0);
  arb_clear(Lu0); arb_clear(dLu0);
  
  arb_clear(aux);
  
  return 0;
}


int compute_big_y_T1(arb_ptr out, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t aux, aux2;
  arb_t term1,term2,term3,term4;
  arb_t factor1,factor2,factor3;
  arb_t const_pi, pm_one;
  arb_t error_tail;
  
  int n;
  
  arb_init(aux); arb_init(aux2);
  arb_init(term1); arb_init(term2); arb_init(term3); arb_init(term4);
  arb_init(factor1); arb_init(factor2); arb_init(factor3);
  arb_init(const_pi); arb_init(pm_one);
  arb_init(error_tail);
  
  arb_const_pi(const_pi,prec);
  arb_zero_pm_one(pm_one);
  
  arb_set_str(aux,"1.5",prec);
  arb_clausen_sin(term1,aux,Ips->x,prec);
  
  arb_mul_si(term2,Ips->x,2,prec);
  arb_clausen_sin(term2,aux,term2,prec);
  arb_div_si(term2,term2,2,prec);
  arb_sqrt_ui(aux,2,prec);
  arb_sub_si(aux,aux,1,prec);
  arb_neg(aux,aux);
  arb_mul(term2,term2,aux,prec);
  
  arb_mul_si(term3,Ips->x,2,prec);
  arb_set_str(aux,"2.5",prec);
  arb_clausen_cos(term3,aux,term3,prec);
  arb_zeta(aux,aux,prec);
  arb_sub(term3,term3,aux,prec);
  arb_sqrt_ui(aux,2,prec);
  arb_sub_si(aux,aux,2,prec);
  arb_div_si(aux,aux,4,prec);
  arb_mul(term3,term3,aux,prec);
  arb_mul(term3,term3,Ips->u+1,prec);
  arb_div(term3,term3,Ips->u,prec);
  arb_neg(term3,term3);
  
  for (n = 1; n <= TANH_MAX; n++){
      arb_set_si(aux,n);
      arb_tanh(factor1,aux,prec);
      arb_sqrt(factor1,factor1,prec);
      arb_sub_si(factor1,factor1,1,prec);
      arb_div_si(factor1,factor1,n,prec);
      arb_sqrt(aux,aux,prec);
      arb_div(factor1,factor1,aux,prec);
      arb_neg(factor1,factor1);
      
      arb_mul_si(aux,Ips->x,n,prec);
      
      arb_cos(factor2,aux,prec);
      arb_sub_si(factor2,factor2,1,prec);
      
      arb_mul(factor2,factor2,Ips->u+1,prec);
      arb_div(factor2,factor2,Ips->u,prec);
      arb_div_si(factor2,factor2,n,prec);
      
      arb_sin(aux,aux,prec);
      arb_add(factor2,factor2,aux,prec);
      
      arb_mul_si(aux,Ips->x,n,prec);
      arb_cos(factor3,aux,prec);
      if (n%2){
        arb_add_si(factor3,factor3,1,prec);
      }
      else {
        arb_sub_si(factor3,factor3,1,prec);
      }
      arb_neg(factor3,factor3);
      
      // multiply and collect the factors
      arb_mul(aux,factor1,factor2,prec);
      arb_mul(aux,aux,factor3,prec);
      
      arb_add(term4,term4,aux,prec);
  }

  // Tail error
  arb_set_str(aux,"-1.5",prec);
  compute_tail_error(error_tail,aux,TANH_MAX+1,prec);
  arb_abs(error_tail,error_tail);
  
  arb_set_str(aux,"-2.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_mul(aux,aux,Ips->u+1,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(error_tail,error_tail,aux,prec);
  
  arb_mul_si(error_tail,error_tail,2,prec);
  arb_mul(error_tail,error_tail,pm_one,prec);
  
  arb_add(out,term1,term2,prec);
  arb_add(out,out,term3,prec);
  arb_add(out,out,term4,prec);
  arb_add(out,out,error_tail,prec);
  arb_div(out,out,Ips->u,prec);
  arb_div(out,out,const_pi,prec);
  
  arb_clear(aux); arb_clear(aux2);
  arb_clear(term1); arb_clear(term2); arb_clear(term3); arb_clear(term4);
  arb_clear(factor1); arb_clear(factor2); arb_clear(factor3);
  arb_clear(const_pi); arb_clear(pm_one);
  arb_clear(error_tail);
  return 1;
}


