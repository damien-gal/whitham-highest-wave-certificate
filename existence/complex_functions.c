/* Blah blah blah. Javi */

#include "complex_functions.h"


int compute_alpha0(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0;
  arb_init(u0);
  
  compute_u0(u0,inp,params,prec);
    
  arb_div_si(out,inp,2,prec);
  arb_div(out,out,u0,prec);
  arb_abs(out,out);
  
  arb_clear(u0);
  
  return 0;
}

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

int compute_alpha2bar(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0,du0;
  arb_init(u0); arb_init(du0);
  
  compute_u0(u0,inp,params,prec);   
  compute_du0(du0,inp,params,prec);
  
  arb_div(out,inp,u0,prec);
  arb_sqr(out,out,prec);
  arb_mul(out,out,du0,prec);
  arb_mul_si(out,out,2,prec);
  arb_abs(out,out);
  
  arb_clear(u0);
  arb_clear(du0);
  
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


void compute_K1(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);

  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IKp(out,argument,prec);

  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  IKp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
}


void compute_K20(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);

  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IKp(out,argument,prec);
  arb_neg(out,out);

  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
  
  IKp(aux,Ips->x,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
}


void compute_K2bar0(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);

  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IIKp(out,argument,prec);

  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IIKp(aux,argument,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  IIKp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  IIKp(aux,Ips->x,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
}


void compute_dyK20(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);

  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(out,argument,prec);
  arb_neg(out,out);

  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_add(out,out,aux,prec);
    
  arb_clear(argument);
  arb_clear(aux);
}


void compute_dyK2bar0(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);

  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IKp(out,argument,prec);

  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  IKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  IKp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
}


void compute_ddyK20(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);

  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(out,argument,prec);
  arb_neg(out,out);

  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
}


void compute_ddyK2bar0(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);

  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(out,argument,prec);

  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
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
  
  // Simple way. Maybe one can optimize later
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


int compute_ddF0(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0,du0,ddu0;
  arb_t Lu0,dLu0,ddLu0;
  arb_t aux,aux2;
  arb_init(u0); arb_init(du0); arb_init(ddu0);
  arb_init(Lu0); arb_init(dLu0); arb_init(ddLu0);
  arb_init(aux); arb_init(aux2);
  
  compute_u0(u0,inp,params,prec); compute_du0(du0,inp,params,prec); compute_ddu0(ddu0,inp,params,prec);    
  compute_Lu0(Lu0,inp,params,prec); compute_dLu0(dLu0,inp,params,prec); compute_ddLu0(ddLu0,inp,params,prec);
  
  
  // First term: 2 derivatives hitting in 1/x, 0 hitting in Lu/u - u
  
  arb_div(aux2,Lu0,u0,prec);
  arb_sub(aux2,aux2,u0,prec);
  arb_div(aux,aux2,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_div(out,aux,inp,prec);
  
  // Second term: 1 derivative hitting in 1/x, 1 hitting in Lu/u - u
  
  arb_mul(aux2,Lu0,du0,prec);
  arb_div(aux2,aux2,u0,prec);
  arb_sub(aux,dLu0,aux2,prec);
  arb_div(aux,aux,u0,prec);
  arb_sub(aux,aux,du0,prec);
  arb_div(aux,aux,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_sub(out,out,aux,prec);
  
  // Third term: 0 derivatives hitting in 1/x, 2 hitting in Lu/u - u
  // The expression of the derivative of Lu/u - u is given by
  // 1/u*(Lu'' + 1/u*(-2*Lu'*u'-Lu*u'' + 2*Lu/u*(u')^2))
  arb_mul(aux2,du0,du0,prec);
  arb_mul(aux2,aux2,Lu0,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_div(aux,aux2,u0,prec);
  arb_mul(aux2,Lu0,ddu0,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_mul(aux2,dLu0,du0,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,u0,prec);
  arb_add(aux,aux,ddLu0,prec);
  arb_div(aux,aux,u0,prec);
  arb_sub(aux,aux,ddu0,prec);
  arb_div_si(aux,aux,2,prec);
  arb_div(aux,aux,inp,prec);
  arb_add(out,out,aux,prec);
  
    
  arb_clear(u0); arb_clear(du0); arb_clear(ddu0);
  arb_clear(Lu0); arb_clear(dLu0); arb_clear(ddLu0);
  
  arb_clear(aux); arb_clear(aux2);
  
  return 0;
}


// TODO(Javi): Fill this up
int compute_dddF0(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t u0,du0,ddu0,dddu0;
  arb_t Lu0,dLu0,ddLu0,dddLu0;
  arb_t aux,aux2;
  arb_init(u0); arb_init(du0); arb_init(ddu0); arb_init(dddu0);
  arb_init(Lu0); arb_init(dLu0); arb_init(ddLu0); arb_init(dddLu0);
  arb_init(aux); arb_init(aux2);
  
  compute_u0(u0,inp,params,prec); compute_du0(du0,inp,params,prec); 
  compute_ddu0(ddu0,inp,params,prec); compute_dddu0(dddu0,inp,params,prec);
  compute_Lu0(Lu0,inp,params,prec); compute_dLu0(dLu0,inp,params,prec); 
  compute_ddLu0(ddLu0,inp,params,prec); compute_dddLu0(dddLu0,inp,params,prec);
  
  
  // First term: 3 derivatives hitting in 1/x, 0 hitting in Lu/u - u
  
  arb_div(aux2,Lu0,u0,prec);
  arb_sub(aux2,aux2,u0,prec);
  arb_div(aux,aux2,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_mul_si(out,aux,-3,prec);
  
  // Second term: 2 derivatives hitting in 1/x, 1 hitting in Lu/u - u
  // The expression of the derivative of Lu/u - u is given by
  // 1/u*(Lu' - Lu*u'/u) - u'
  
  arb_mul(aux2,Lu0,du0,prec);
  arb_div(aux2,aux2,u0,prec);
  arb_sub(aux,dLu0,aux2,prec);
  arb_div(aux,aux,u0,prec);
  arb_sub(aux,aux,du0,prec);
  arb_div(aux,aux,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_mul_si(aux,aux,3,prec);
  arb_add(out,out,aux,prec);
  
  // Third term: 1 derivative hitting in 1/x, 2 hitting in Lu/u - u
  // The expression of the derivative of Lu/u - u is given by
  // 1/u*(Lu'' + 1/u*(-2*Lu'*u'-Lu*u'' + 2*Lu/u*(u')^2)) - u''
  arb_mul(aux2,du0,du0,prec);
  arb_mul(aux2,aux2,Lu0,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_div(aux,aux2,u0,prec);
  arb_mul(aux2,Lu0,ddu0,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_mul(aux2,dLu0,du0,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,u0,prec);
  arb_add(aux,aux,ddLu0,prec);
  arb_div(aux,aux,u0,prec);
  arb_sub(aux,aux,ddu0,prec);
  arb_div_si(aux,aux,2,prec);
  arb_mul_si(aux,aux,-3,prec);
  arb_div(aux,aux,inp,prec);
  arb_div(aux,aux,inp,prec);
  arb_add(out,out,aux,prec);
  
  
  // Fourth term: 0 derivatives hitting in 1/x, 3 hitting in Lu/u - u
  // The expression of the derivative of Lu/u - u is given by
  // 1/u*{Lu''' + 1/u*[-3*u'*Lu''-3*Lu'*u''-u'''*Lu
  //      +1/u*(6*(u')^2*Lu' + 6*u''*u'*Lu + 1/u*{-6*(u')^3*Lu})]} - u'''
  arb_mul(aux2,du0,du0,prec);
  arb_mul(aux2,aux2,du0,prec);
  arb_mul(aux2,aux2,Lu0,prec);
  arb_mul_si(aux2,aux2,-6,prec);
  arb_div(aux,aux2,u0,prec);
  
  arb_mul(aux2,ddu0,du0,prec);
  arb_mul(aux2,aux2,Lu0,prec);
  arb_mul_si(aux2,aux2,6,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_mul(aux2,du0,du0,prec);
  arb_mul(aux2,aux2,dLu0,prec);
  arb_mul_si(aux2,aux2,6,prec);
  arb_add(aux,aux,aux2,prec);
  arb_div(aux,aux,u0,prec);
  
  arb_mul(aux2,dddu0,Lu0,prec);
  arb_sub(aux,aux,aux2,prec);
  
  arb_mul(aux2,dLu0,ddu0,prec);
  arb_mul_si(aux2,aux2,-3,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_mul(aux2,ddLu0,du0,prec);
  arb_mul_si(aux2,aux2,-3,prec);
  arb_add(aux,aux,aux2,prec);
  arb_div(aux,aux,u0,prec);
  
  arb_add(aux,aux,dddLu0,prec);
  arb_div(aux,aux,u0,prec);
  
  arb_sub(aux,aux,dddu0,prec);
  
  arb_div(aux,aux,inp,prec);
  arb_div_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
    
  arb_clear(u0); arb_clear(du0); arb_clear(ddu0); arb_clear(dddu0);
  arb_clear(Lu0); arb_clear(dLu0); arb_clear(ddLu0); arb_clear(dddLu0);
  
  arb_clear(aux); arb_clear(aux2);
  
  return 0;
}

int compute_F1(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t F0, dF0;
  arb_init(F0); arb_init(dF0);
  
  compute_F0(F0,inp,params,prec);
  compute_dF0(dF0,inp,params,prec);
  
  arb_mul(out,dF0,inp,prec);
  arb_add(out,out,F0,prec);
  
  arb_clear(F0); arb_clear(dF0);
  
  return 0;
}

int compute_dF1(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t dF0, ddF0;
  arb_init(dF0); arb_init(ddF0);

  compute_dF0(dF0,inp,params,prec);
  compute_ddF0(ddF0,inp,params,prec);
  
  arb_mul(out,ddF0,inp,prec);
  arb_mul_si(dF0,dF0,2,prec);
  arb_add(out,out,dF0,prec);
  
  arb_clear(dF0); arb_clear(ddF0);
  
  return 0;
}

int compute_ddF1(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t ddF0, dddF0;
  arb_init(ddF0); arb_init(dddF0);

  compute_ddF0(ddF0,inp,params,prec);
  compute_dddF0(dddF0,inp,params,prec);
  
  arb_mul(out,dddF0,inp,prec);
  arb_mul_si(ddF0,ddF0,3,prec);
  arb_add(out,out,ddF0,prec);
  
  arb_clear(ddF0); arb_clear(dddF0);
  
  return 0;
}


int compute_F2(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t F1, dF1;
  arb_init(F1); arb_init(dF1);

  compute_F1(F1,inp,params,prec);
  compute_dF1(dF1,inp,params,prec);
  
  arb_mul(out,dF1,inp,prec);
  arb_add(out,out,F1,prec);
  
  arb_clear(F1); arb_clear(dF1);
  
  return 0;
}

int compute_dF2(arb_ptr out, const arb_t inp, void* params, slong prec){
  arb_t dF1, ddF1;
  arb_init(dF1); arb_init(ddF1);

  compute_dF1(dF1,inp,params,prec);
  compute_ddF1(ddF1,inp,params,prec);
  
  arb_mul(out,ddF1,inp,prec);
  arb_add(out,out,dF1,prec);
  
  arb_clear(dF1); arb_clear(ddF1);
  
  return 0;
}


// TODO(Javi): untested for tail error
int compute_big_y_T0(arb_ptr out, void* params, slong prec){
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
  arb_neg(term2,term2);
  
  arb_mul_si(term3,Ips->x,2,prec);
  arb_set_str(aux,"2.5",prec);
  arb_clausen_cos(term3,aux,term3,prec);
  arb_zeta(aux,aux,prec);
  arb_sub(term3,term3,aux,prec);
  arb_sqrt_ui(aux,2,prec);
  arb_sub_si(aux,aux,2,prec);
  arb_div_si(aux,aux,4,prec);
  arb_div(aux,aux,Ips->x,prec);
  arb_mul(term3,term3,aux,prec);
    
  for (n = 1; n <= TANH_MAX; n++){
  
      arb_set_si(aux,n);
      arb_tanh(factor1,aux,prec);
      arb_sqrt(factor1,factor1,prec);
      arb_sub_si(factor1,factor1,1,prec);
      arb_div_si(factor1,factor1,n*n,prec);
      arb_sqrt(aux,aux,prec);
      arb_div(factor1,factor1,aux,prec);
      
      arb_mul_si(aux,Ips->x,n,prec);
      
      arb_cos(factor2,aux,prec);
      arb_sub_si(factor2,factor2,1,prec);
      
      arb_cos(factor3,aux,prec);
      arb_neg(factor3,factor3);
      
      arb_sin(aux2,aux,prec);
      arb_mul(aux2,aux2,aux,prec);
      arb_sub(factor3,factor3,aux2,prec);
      if (n%2){
        arb_sub_si(factor3,factor3,1,prec);
      }
      else {
        arb_add_si(factor3,factor3,1,prec);
      }
      
      // multiply and collect the factors
      arb_mul(aux,factor1,factor2,prec);
      arb_mul(aux,aux,factor3,prec);
      arb_div(aux,aux,Ips->x,prec);
      
      arb_add(term4,term4,aux,prec);
  }
  
  // Tail error
  arb_set_str(aux,"-1.5",prec);
  compute_tail_error(error_tail,aux,TANH_MAX+1,prec);
  arb_abs(error_tail,error_tail);
  
  arb_set_str(aux,"-2.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_div(aux,aux,Ips->x,prec);
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


int compute_big_y_T2(arb_ptr out, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t aux, aux2, aux3;
  arb_t const_pi, pm_one;
  
  arb_t SE;
  arb_t error_SE;
  
  arb_t I0, I1, I2;
  
  int n;
  int M = 50; // TODO(Javi): CHECK THIS WELL!!!!
  
  arb_init(aux); arb_init(aux2); arb_init(aux3);
  arb_init(const_pi); arb_init(pm_one);
  arb_init(SE);
  arb_init(error_SE);
  arb_init(I0); arb_init(I1); arb_init(I2);
  
  arb_const_pi(const_pi,prec);
  arb_zero_pm_one(pm_one);

  arb_set_str(out,"0",prec);
  arb_set_str(SE,"0",prec);
  arb_set_str(I0,"0",prec);
  arb_set_str(I1,"0",prec);
  arb_set_str(I2,"0",prec);  

  // SE
  for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux,Ips->x,n,prec);
      arb_cos(aux3,aux,prec);
      
      arb_sin(aux2,aux,prec);
      arb_mul(aux2,aux2,Ips->u+1,prec);
      arb_div(aux2,aux2,Ips->u,prec);
      arb_mul_si(aux2,aux2,2,prec);
      arb_div_si(aux2,aux2,n,prec);
      arb_sub(aux3,aux3,aux2,prec);
      
      arb_sqr(aux,Ips->u+1,prec);
      arb_mul_si(aux,aux,2,prec);
      arb_mul(aux2,Ips->u,Ips->u+2,prec);
      arb_sub(aux,aux,aux2,prec);
      
      arb_mul_si(aux2,Ips->x,n,prec);
      arb_cos(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      arb_neg(aux2,aux2);
      arb_mul(aux,aux,aux2,prec);
      
      arb_div(aux,aux,Ips->u,prec);
      arb_div(aux,aux,Ips->u,prec);
      arb_div_si(aux,aux,n*n,prec);
      arb_add(aux3,aux3,aux,prec);
   
            
      arb_mul_si(aux,const_pi,n,prec);
      arb_cosintegral(aux2,aux,prec);
      arb_mul_si(aux,Ips->x,n,prec);
      arb_cosintegral(aux,aux,prec);
      arb_sub(aux2,aux2,aux,prec);
      arb_mul(aux3,aux3,aux2,prec);
      
      arb_set_si(aux2,n);
      arb_tanh(aux2,aux2,prec);
      arb_sqrt(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      arb_set_si(aux,n);
      arb_sqrt(aux,aux,prec);
      arb_div(aux2,aux2,aux,prec);
      arb_mul(aux3,aux3,aux2,prec);
      
      arb_add(SE,SE,aux3,prec);
  }
  arb_mul(SE,SE,Ips->x,prec);
  arb_div(SE,SE,const_pi,prec);
  arb_div(SE,SE,Ips->u,prec);
  
  // Error-SE
  arb_sqr(aux,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,Ips->u,Ips->u+2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_set_str(aux2,"-2.5",prec);
  compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_abs(aux,aux);
  
  arb_set_str(aux2,"-1.5",prec);
  compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_mul(aux2,aux2,Ips->u+1,prec);
  arb_div(aux2,aux2,Ips->u,prec);
  arb_abs(aux2,aux2);
  arb_add(aux,aux,aux2,prec);
  
  arb_set_str(aux2,"-0.5",prec);
  compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_mul(error_SE,aux,Ips->x,prec);
  arb_div(error_SE,error_SE,Ips->u,prec);
  
  arb_div(aux,const_pi,Ips->x,prec);
  arb_log(aux,aux,prec);
  arb_mul(error_SE,error_SE,aux,prec);
  arb_div(error_SE,error_SE,const_pi,prec);
  
  arb_mul(error_SE,error_SE,pm_one,prec);
  

  // I0
  for (n = 1; n <= M; n++){ 
      arb_mul_si(aux,Ips->x,n,prec);
      arb_cos(aux3,aux,prec);
      
      arb_sin(aux2,aux,prec);
      arb_mul(aux2,aux2,Ips->u+1,prec);
      arb_div(aux2,aux2,Ips->u,prec);
      arb_mul_si(aux2,aux2,2,prec);
      arb_div_si(aux2,aux2,n,prec);
      arb_sub(aux3,aux3,aux2,prec);
      
      arb_sqr(aux,Ips->u+1,prec);
      arb_mul_si(aux,aux,2,prec);
      arb_mul(aux2,Ips->u,Ips->u+2,prec);
      arb_sub(aux,aux,aux2,prec);
      
      arb_mul_si(aux2,Ips->x,n,prec);
      arb_cos(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      arb_neg(aux2,aux2);
      arb_mul(aux,aux,aux2,prec);
      
      arb_div(aux,aux,Ips->u,prec);
      arb_div(aux,aux,Ips->u,prec);
      arb_div_si(aux,aux,n*n,prec);
      arb_add(aux3,aux3,aux,prec);
   
            
      arb_mul_si(aux,const_pi,n,prec);
      arb_cosintegral(aux2,aux,prec);
      arb_mul_si(aux,Ips->x,n,prec);
      arb_cosintegral(aux,aux,prec);
      arb_sub(aux2,aux2,aux,prec);

      arb_mul_si(aux,Ips->x,n,prec);
      arb_cos(aux,aux,prec);
      arb_sub_si(aux,aux,1,prec);
      arb_div_si(aux,aux,n*n,prec);
      arb_div(aux,aux,Ips->x,prec);
      arb_div(aux,aux,Ips->x,prec);
      arb_sub(aux2,aux2,aux,prec);
      
      arb_mul_si(aux,Ips->x,n,prec);
      arb_sin(aux,aux,prec);
      arb_div_si(aux,aux,n,prec);
      arb_div(aux,aux,Ips->x,prec);
      arb_add(aux2,aux2,aux,prec);
      
      if (n%2){ // -2/(n^2*pi^2)
        arb_set_si(aux,2);
        arb_div_si(aux,aux,n*n,prec);
        arb_div(aux,aux,const_pi,prec);
        arb_div(aux,aux,const_pi,prec);
        arb_sub(aux2,aux2,aux,prec);
      }
      
      arb_mul(aux3,aux3,aux2,prec);
      
      arb_set_si(aux2,n);
      arb_sqrt(aux2,aux2,prec);
      arb_div(aux3,aux3,aux2,prec);
      arb_add(I0,I0,aux3,prec);
  }
  arb_mul(I0,I0,Ips->x,prec);
  arb_div(I0,I0,const_pi,prec);
  arb_div(I0,I0,Ips->u,prec);
  
  // TODO(Javi): TEST THIS CAREFULLY!
  
  // 3rd and 4th lines:
  // 3rd line
  arb_sqr(aux,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,Ips->u,Ips->u+2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_div(I1,aux,Ips->u,prec);
    
  arb_set_str(aux,"3.5",prec);
  arb_mul_si(aux2,Ips->x,2,prec);
  arb_clausen_sin(aux,aux,aux2,prec);
  arb_mul(aux,aux,Ips->x,prec);
  arb_mul(aux,aux,const_pi,prec);
  arb_mul(aux,aux,const_pi,prec);
  printf("aux = %s\n",arb_get_str(aux,prec,1));

  arb_set_str(aux2,"3.5",prec);
  arb_clausen_sin(aux2,aux2,Ips->x,prec);
  arb_mul(aux2,aux2,Ips->x,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul_si(aux2,aux2,2,prec);
  printf("aux2 = %s\n",arb_get_str(aux2,prec,1));
  arb_sub(aux,aux,aux2,prec);
  
  arb_sqr(aux2,const_pi,prec);
  arb_sqr(aux3,Ips->x,prec);
  arb_sub(aux2,aux2,aux3,prec);
  arb_mul_si(aux2,aux2,4,prec);
  arb_set_str(aux3,"4.5",prec);
  arb_clausen_cos(aux3,aux3,Ips->x,prec);
  arb_mul(aux2,aux2,aux3,prec);
  printf("aux2 = %s\n",arb_get_str(aux2,prec,1));
  arb_add(aux,aux,aux2,prec);
    
  // 4th line
  arb_set_str(aux3,"4.5",prec);
  arb_mul_si(aux2,Ips->x,2,prec);
  arb_clausen_cos(aux2,aux3,aux2,prec);
  arb_mul(aux2,aux2,Ips->x,prec);
  arb_mul(aux2,aux2,Ips->x,prec);
  arb_set_si(aux3,32);
  arb_sqrt(aux3,aux3,prec);
  arb_div(aux2,aux2,aux3,prec);
  printf("aux2 = %s\n",arb_get_str(aux2,prec,1));
  arb_add(aux,aux,aux2,prec);
  
  arb_set_str(aux3,"4.5",prec);
  arb_mul_si(aux2,Ips->x,2,prec);
  arb_clausen_cos(aux2,aux3,aux2,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  printf("aux2 = %s\n",arb_get_str(aux2,prec,1));
  arb_sub(aux,aux,aux2,prec);
    
  arb_set_si(aux3,32);
  arb_sqrt(aux3,aux3,prec);
  arb_inv(aux3,aux3,prec);
  arb_sub_si(aux3,aux3,4,prec);
  arb_mul(aux3,aux3,Ips->x,prec);
  arb_mul(aux3,aux3,Ips->x,prec);
  arb_sqr(aux2,const_pi,prec);
  arb_mul_si(aux2,aux2,3,prec);
  arb_add(aux3,aux2,aux3,prec);
  arb_set_str(aux2,"4.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux3,aux3,aux2,prec);
  printf("aux3 = %s\n",arb_get_str(aux3,prec,1));
  arb_sub(aux,aux,aux3,prec);
  printf("aux = %s\n",arb_get_str(aux,prec,1));
  
  arb_mul(I1,I1,aux,prec);
  
  // 2nd line
  arb_set_str(aux,"2.5",prec);
  arb_mul_si(aux2,Ips->x,2,prec);
  arb_clausen_cos(aux,aux,aux2,prec);
  arb_mul(aux,aux,Ips->x,prec);
  arb_mul(aux,aux,const_pi,prec);
  arb_mul(aux,aux,const_pi,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_neg(aux,aux);
  
  arb_sqr(aux2,const_pi,prec);
  arb_sqr(aux3,Ips->x,prec);
  arb_mul_si(aux3,aux3,2,prec);
  arb_sub(aux2,aux2,aux3,prec);
  arb_mul_si(aux2,aux2,4,prec);
  arb_set_str(aux3,"3.5",prec);
  arb_clausen_sin(aux3,aux3,Ips->x,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_set_str(aux3,"3.5",prec);
  arb_mul_si(aux2,Ips->x,2,prec);
  arb_clausen_sin(aux2,aux3,aux2,prec);
  arb_mul(aux2,aux2,Ips->x,prec);
  arb_mul(aux2,aux2,Ips->x,prec);
  arb_set_si(aux3,2);
  arb_sqrt(aux3,aux3,prec);
  arb_div(aux2,aux2,aux3,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_set_str(aux3,"3.5",prec);
  arb_mul_si(aux2,Ips->x,2,prec);
  arb_clausen_sin(aux2,aux3,aux2,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_sub(aux,aux,aux2,prec);
  
  arb_set_str(aux2,"2.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux2,aux2,Ips->x,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_mul(aux,aux,Ips->u+1,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_add(I1,I1,aux,prec);
  
  // 1st line
  arb_sqr(aux2,const_pi,prec);
  arb_sqr(aux3,Ips->x,prec);
  arb_mul_si(aux3,aux3,2,prec);
  arb_sub(aux2,aux2,aux3,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_set_str(aux3,"2.5",prec);
  arb_clausen_cos(aux3,aux3,Ips->x,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_sub(I1,I1,aux2,prec);
  
  arb_sqr(aux2,const_pi,prec);
  arb_sqr(aux3,Ips->x,prec);
  arb_set_si(aux,2);
  arb_sqrt(aux,aux,prec);
  arb_div(aux3,aux3,aux,prec);
  arb_sub(aux2,aux2,aux3,prec);
  arb_set_str(aux3,"2.5",prec);
  arb_mul_si(aux,Ips->x,2,prec);
  arb_clausen_cos(aux3,aux3,aux,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_add(I1,I1,aux2,prec);
  
  arb_set_str(aux3,"1.5",prec);
  arb_mul_si(aux2,Ips->x,2,prec);
  arb_clausen_sin(aux3,aux3,aux2,prec);
  arb_mul(aux3,aux3,Ips->x,prec);
  arb_mul(aux3,aux3,const_pi,prec);
  arb_mul(aux3,aux3,const_pi,prec);
  arb_sub(I1,I1,aux3,prec);

  arb_set_str(aux2,"2.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_mul(aux2,aux2,const_pi,prec);
  arb_add(I1,I1,aux2,prec);
  
  arb_div_si(I1,I1,2,prec);
  arb_div(I1,I1,Ips->x,prec);
  arb_div(I1,I1,const_pi,prec);
  arb_div(I1,I1,const_pi,prec);
  arb_div(I1,I1,const_pi,prec);
  arb_div(I1,I1,Ips->u,prec);
  
  
  // Calculation of I2
  arb_sqr(aux,const_pi,prec);
  arb_inv(aux,aux,prec);
  arb_sqr(aux2,Ips->x,prec);
  arb_inv(aux2,aux2,prec);
  arb_sub(I2,aux2,aux,prec);
    
  arb_sqr(aux,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,Ips->u,Ips->u+2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_div_si(aux,aux,7,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_set_str(aux2,"-3.5",prec);
  arb_set_si(aux3,M);
  arb_pow(aux2,aux3,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  
  arb_set_str(aux2,"-2.5",prec);
  arb_set_si(aux3,M);
  arb_pow(aux2,aux3,aux2,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_div_si(aux2,aux2,5,prec);
  arb_mul(aux2,aux2,Ips->u+1,prec);
  arb_div(aux2,aux2,Ips->u,prec);
  arb_abs(aux2,aux2);
  arb_add(aux,aux,aux2,prec);
  
  arb_set_str(aux2,"-1.5",prec);
  arb_set_si(aux3,M);
  arb_pow(aux2,aux3,aux2,prec);
  arb_div_si(aux2,aux2,3,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_mul(I2,I2,aux,prec);
  
  arb_mul_si(I2,I2,4,prec);
  arb_mul(I2,I2,Ips->x,prec);
  arb_div(I2,I2,Ips->u,prec);
  arb_div(I2,I2,const_pi,prec);
  arb_mul(I2,I2,pm_one,prec);
  
  arb_add(out,SE,error_SE,prec);
  arb_add(out,out,I0,prec);
  arb_add(out,out,I1,prec);
  arb_add(out,out,I2,prec);
  
  printf("SE = %s\n",arb_get_str(SE,prec,1));
  printf("error_SE = %s\n",arb_get_str(error_SE,prec,1));
  printf("I0 = %s\n",arb_get_str(I0,prec,1));
  printf("I1 = %s\n",arb_get_str(I1,prec,1));
  printf("I2 = %s\n",arb_get_str(I2,prec,1));
  
  arb_clear(I0); arb_clear(I1); arb_clear(I2);
  arb_clear(SE); arb_clear(error_SE);
  
  arb_clear(aux); arb_clear(aux2); arb_clear(aux3);
  arb_clear(const_pi);
  arb_clear(pm_one);
  return 1;
}


int integrand_T0(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);
  
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(out,argument,prec);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(out,out,inp,prec);
  
  arb_div_si(out,out,2,prec);
  arb_mul(out,out,Ips->x,prec);
  arb_div(out,out,Ips->u,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 0;
}

// TODO(Javi): Do I need this one???!?
/*int integrand_dT0(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);
  
  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(out,argument,prec);
  
  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(out,out,inp,prec);
  arb_mul(out,out,Ips->x,prec);
  
  
  
  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  
  arb_clear(argument);
  arb_clear(aux);
  return 0;
}*/


int integrand_ddT0(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);
  
    
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(out,argument,prec);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul_si(out,out,2,prec);
  
  
  // TODO(Javi): Collapse these 3 into 1 and multiply by Ips->x at the end, not term by term
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  ddKp(aux,argument,prec);
  arb_mul(aux,aux,Ips->x,prec);
  arb_mul(aux,aux,inp,prec);
  arb_add(out,out,aux,prec);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  ddKp(aux,argument,prec);
  arb_mul(aux,aux,Ips->x,prec);
  arb_mul(aux,aux,inp,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(argument,inp,Ips->x,prec);
  ddKp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux,aux,Ips->x,prec);
  arb_mul(aux,aux,inp,prec);
  arb_sub(out,out,aux,prec);
  
  
  arb_mul(out,out,Ips->x,prec);
  
  
  arb_div_si(out,out,2,prec);
  arb_mul(out,out,Ips->x,prec);
  arb_div(out,out,Ips->u,prec);

  arb_clear(argument);
  arb_clear(aux);
  return 0;
}


int integrand_T1(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(out,argument,prec);
  
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
  
  compute_K1(aux,inp,params,prec);  
  arb_mul(aux,aux,Ips->u+1,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_sub(out,out,aux,prec);
  
  arb_div_si(out,out,2,prec);
  arb_mul(out,out,Ips->x,prec);
  arb_div(out,out,Ips->u,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 0;
}

// TODO(Javi): Change the "+1" and "-1" as in T0
int integrand_ddT1(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;

  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);
  
  
  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(out,argument,prec);
  
  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);
    
  arb_mul(argument,inp,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_mul(out,out,Ips->u+1,prec);
  arb_div(out,out,Ips->u,prec);
  
  
  arb_set_si(argument,1);
  arb_sub(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  ddKp(aux,argument,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_si(argument,1);
  arb_add(argument,argument,inp,prec);
  arb_mul(argument,argument,Ips->x,prec);
  ddKp(aux,argument,prec);
  arb_sub(out,out,aux,prec);  
  
  
  arb_mul(out,out,Ips->x,prec);
  arb_mul(out,out,Ips->x,prec);
  
  
  arb_div_si(out,out,2,prec);
  arb_mul(out,out,Ips->x,prec);
  arb_div(out,out,Ips->u,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 0;
}


int integrand_T2(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux, aux2;
  
  arb_init(argument);
  arb_init(aux); arb_init(aux2);
 
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(out,argument,prec);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_add(out,out,aux,prec);
  
  Kp(aux,Ips->x,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  
   
  arb_div(aux,Ips->u+1,Ips->u,prec);
  compute_K20(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  
  arb_mul(aux,Ips->u+1,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,Ips->u,Ips->u+2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_div(aux,aux,Ips->u,prec);
  compute_K2bar0(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  
  arb_div(out,out,inp,prec);
  
  arb_div_si(out,out,2,prec);
  arb_mul(out,out,Ips->x,prec);
  arb_div(out,out,Ips->u,prec);
  
  
  arb_clear(argument);
  arb_clear(aux); arb_clear(aux2);
  return 0;
}

// TODO(Javi): untested.
// VERY VERY INEFFICIENT
int integrand_ddT2(arb_ptr out, const arb_t inp, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t argument;
  arb_t aux, aux2, aux3;
  
  arb_init(argument);
  arb_init(aux); arb_init(aux2); arb_init(aux3);
 
  // Two derivatives in 1/w, 0 in the big parenthesis 
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux3,argument,prec);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  Kp(aux,argument,prec);
  arb_add(aux3,aux3,aux,prec);
  
  Kp(aux,Ips->x,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(aux3,aux3,aux,prec);
  
   
  arb_div(aux,Ips->u+1,Ips->u,prec);
  compute_K20(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(aux3,aux3,aux,prec);
  
  
  arb_mul(aux,Ips->u+1,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,Ips->u,Ips->u+2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_div(aux,aux,Ips->u,prec);
  compute_K2bar0(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(aux3,aux3,aux,prec);
  
  arb_div(aux3,aux3,inp,prec);
  arb_div(aux3,aux3,inp,prec);
  arb_div(aux3,aux3,inp,prec);
  arb_mul_si(out,aux3,2,prec);
  

  // 2*one derivative in 1/w, one in the big parenthesis
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(aux3,argument,prec);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  dKp(aux,argument,prec);
  arb_sub(aux3,aux3,aux,prec);
   
  arb_div(aux,Ips->u+1,Ips->u,prec);
  compute_dyK20(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(aux3,aux3,aux,prec);
  
  
  arb_mul(aux,Ips->u+1,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,Ips->u,Ips->u+2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_div(aux,aux,Ips->u,prec);
  compute_dyK2bar0(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(aux3,aux3,aux,prec);
  
  
  arb_div(aux3,aux3,inp,prec);
  arb_div(aux3,aux3,inp,prec);
  arb_mul_si(aux3,aux3,-2,prec);
  arb_mul(aux3,aux3,Ips->x,prec);
  arb_add(out,out,aux3,prec);
  
  
  // zero derivatives in 1/w, two in the big parenthesis
  arb_add_si(argument,inp,1,prec);
  arb_mul(argument,argument,Ips->x,prec);
  ddKp(aux3,argument,prec);
  
  arb_sub_si(argument,inp,1,prec);
  arb_neg(argument,argument);
  arb_mul(argument,argument,Ips->x,prec);
  ddKp(aux,argument,prec);
  arb_add(aux3,aux3,aux,prec);
     
  arb_div(aux,Ips->u+1,Ips->u,prec);
  compute_ddyK20(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(aux3,aux3,aux,prec);
  
  
  arb_mul(aux,Ips->u+1,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_mul(aux2,Ips->u,Ips->u+2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_div(aux,aux,Ips->u,prec);
  compute_ddyK2bar0(aux2,inp,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(aux3,aux3,aux,prec);
  
  
  arb_div(aux3,aux3,inp,prec);
  arb_mul(aux3,aux3,Ips->x,prec);
  arb_mul(aux3,aux3,Ips->x,prec);
  arb_add(out,out,aux3,prec);
  
  
  
  arb_div_si(out,out,2,prec);
  arb_mul(out,out,Ips->x,prec);
  arb_div(out,out,Ips->u,prec);
  
  
  arb_clear(argument);
  arb_clear(aux); arb_clear(aux2); arb_clear(aux3);
  return 0;
}

int singularity_T0_0(arb_ptr out, const arb_t delta, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t aux, aux2;
  arb_t pi;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // Second line
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m*m*m);
    arb_sqrt(aux2,aux2,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  arb_set_str(aux,"1.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_add(out,out,aux,prec);
  
  arb_div(out,out,pi,prec);
  arb_div_si(out,out,2,prec);
  arb_sqrt(aux,delta,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(aux,"2.5",prec);
  arb_pow(aux,Ips->x,aux,prec);
  arb_mul(out,out,aux,prec);
  
  // First line
  
  arb_set_str(aux,"2.5",prec);
  arb_zeta(aux,aux,prec);
  arb_div_si(aux,aux,2,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_set_str(aux2,"2",prec);
  arb_div(aux2,aux2,pi,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  compute_C_Delta(aux2,delta,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_sqrt(aux2,delta,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_str(aux2,"2.5",prec);
  arb_pow(aux2,Ips->x,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  compute_c_delta(aux,delta,prec);
  arb_sqrt(aux2,delta,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,3,prec);
  arb_add_si(aux,aux,8,prec);
  arb_mul_si(aux2,pi,2,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_div(aux,aux,aux2,prec);
  arb_div_si(aux,aux,6,prec);
  arb_add(out,out,aux,prec);
  
  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(aux,"1.5",prec);
  arb_pow(aux,delta,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,Ips->u,prec);
    
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  return 1;
}


int singularity_T0_1(arb_ptr out, const arb_t delta, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t aux, aux2;
  arb_t pi;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // Second line
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m);
    arb_sqrt(aux2,aux2,prec);
    
    arb_div(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  arb_set_str(aux,"-0.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_add(out,out,aux,prec);
  
  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_mul_si(out,out,3,prec);
  
  arb_set_str(aux,"2.5",prec);
  arb_zeta(aux,aux,prec);
  arb_mul(aux,aux,delta,prec);
  arb_set_str(aux2,"2.5",prec);
  arb_pow(aux2,Ips->x,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_div_si(aux,aux,3,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_add(out,out,aux,prec);
  
  // First line
  
  arb_add_si(out,out,2,prec);
  arb_set_si(aux,2);
  arb_sqrt(aux,aux,prec);
  arb_sub(out,out,aux,prec);
  arb_set_si(aux,2);
  arb_sub(aux,aux,delta,prec);
  arb_sqrt(aux,aux,prec);
  arb_add(out,out,aux,prec);
  arb_set_si(aux,1);
  arb_sub(aux,aux,delta,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  arb_sqrt(aux,delta,prec);
  arb_mul_si(aux,aux,3,prec);
  arb_add(out,out,aux,prec);
    
  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_mul(out,out,delta,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,Ips->u,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
    
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  return 1;
}


int singularity_T1_1(arb_ptr out, const arb_t delta, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t aux, aux2, aux3;
  arb_t pi;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // Second line
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m*m*m);
    arb_sqrt(aux2,aux2,prec);
    
    arb_div(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  arb_set_str(aux,"-1.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"1.5",prec);
  arb_zeta(aux,aux,prec);
  arb_add(out,out,aux,prec);
  
  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_mul(out,out,Ips->u+1,prec);
  arb_div(out,out,Ips->u,prec);
  arb_mul_si(out,out,2,prec);
  
  // First line
  
  arb_set_str(aux,"0",prec);
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux2,aux2,1,prec);
    
    arb_set_si(aux3,m);
    arb_sqrt(aux3,aux3,prec);
    
    arb_div(aux2,aux2,aux3,prec);
    arb_sub(aux,aux,aux2,prec);
  }
  arb_set_str(aux2,"-0.5",prec);
  compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_sqrt(aux2,Ips->x,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_si(aux2,2);
  arb_div(aux2,aux2,pi,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"2",prec);
  arb_sub(aux,aux,delta,prec);
  arb_sqrt(aux,aux,prec);
  arb_add(out,out,aux,prec);
  arb_sqrt(aux,delta,prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"2",prec);
  arb_sqrt(aux,aux,prec);
  arb_sub(out,out,aux,prec);
  
  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_mul(out,out,delta,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,Ips->u,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
    
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  return 1;
}



int singularity_T2_0(arb_ptr out, const arb_t delta, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t aux, aux2, aux3;
  arb_t pi;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // Third line
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m);
    arb_sqrt(aux2,aux2,prec);
    
    arb_div(aux,aux,aux2,prec);
    
    arb_mul(aux2,Ips->u+1,Ips->u+1,prec);
    arb_mul_si(aux2,aux2,2,prec);
    arb_div(aux2,aux2,Ips->u,prec);
    arb_sub(aux2,aux2,Ips->u+2,prec);
    arb_abs(aux2,aux2);
    arb_div(aux2,aux2,Ips->u,prec);
    
    arb_div(aux3,Ips->u+1,Ips->u,prec);
    arb_mul_si(aux3,aux3,m,prec);
    arb_add(aux2,aux2,aux3,prec);
    
    arb_set_si(aux3,m*m);
    arb_div_si(aux3,aux3,2,prec);
    arb_add(aux2,aux2,aux3,prec);
    
    arb_mul(aux,aux,aux2,prec);
    
    arb_sub(out,out,aux,prec);
  }
  arb_set_str(aux,"1.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_div_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"0.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_mul(aux,aux,Ips->u+1,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"-0.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_mul(aux2,Ips->u+1,Ips->u+1,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_div(aux2,aux2,Ips->u,prec);
  arb_sub(aux2,aux2,Ips->u+2,prec);
  arb_abs(aux2,aux2);
  arb_div(aux2,aux2,Ips->u,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);

  arb_set_str(aux,"2.5",prec);
  arb_pow(aux,Ips->x,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,2,prec);

  // Second line
  arb_set_str(aux,"0.5",prec);
  arb_zeta(aux,aux,prec);
  arb_abs(aux,aux);
  arb_mul(aux2,Ips->u+1,Ips->u+1,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_div(aux2,aux2,Ips->u,prec);
  arb_sub(aux2,aux2,Ips->u+2,prec);
  arb_abs(aux2,aux2);
  arb_div(aux2,aux2,Ips->u,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_set_si(aux2,2);
  arb_div(aux2,aux2,pi,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_str(aux2,"3.5",prec);
  arb_pow(aux2,Ips->x,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  // First line
  arb_set_str(aux,"-1.5",prec);
  arb_zeta(aux,aux,prec);
  arb_abs(aux,aux);
  arb_set_si(aux2,2);
  arb_div(aux2,aux2,pi,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_str(aux2,"2.5",prec);
  arb_pow(aux2,Ips->x,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_div(aux,Ips->u+1,Ips->u,prec);
  arb_mul(aux,aux,pi,prec);
  arb_mul(aux,aux,pi,prec);
  arb_mul_si(aux,aux,5,prec);
  arb_add_si(aux,aux,1,prec);
  arb_set_str(aux2,"2.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_set_str(aux2,"2.5",prec);
  arb_pow(aux2,Ips->x,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"0.25",prec);
  arb_add(out,out,aux,prec);
  
  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_div(out,out,Ips->u,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_mul(out,out,delta,prec);
  arb_mul(out,out,delta,prec);
  arb_div_si(out,out,2,prec);
  
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  return 1;
}



int singularity_T2_1(arb_ptr out, const arb_t delta, void* params, slong prec){
  Integration_params_struct* Ips = params;
  
  arb_t aux, aux2, aux3;
  arb_t pi;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // Third line
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m*m*m);
    arb_sqrt(aux2,aux2,prec);
    
    arb_div(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  arb_set_str(aux,"-1.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"1.5",prec);
  arb_zeta(aux,aux,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(aux,Ips->u+1,Ips->u+1,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_div(aux,aux,Ips->u,prec);
  arb_sub(aux,aux,Ips->u+2,prec);
  arb_abs(aux,aux);
  arb_add(aux,aux,Ips->u+1,prec);
  arb_mul(out,out,aux,prec);

  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div(out,out,Ips->u,prec);
  arb_mul_si(out,out,2,prec);
  
  // Second line
  
  arb_set_str(aux,"0",prec);
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux2,aux2,1,prec);
    
    arb_set_si(aux3,m);
    arb_sqrt(aux3,aux3,prec);
    
    arb_div(aux2,aux2,aux3,prec);
    arb_sub(aux,aux,aux2,prec);
  }
  arb_set_str(aux2,"-0.5",prec);
  compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_sqrt(aux2,Ips->x,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_si(aux2,2);
  arb_div(aux2,aux2,pi,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"2.5",prec);
  arb_zeta(aux,aux,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_mul(aux2,delta,delta,prec);
  arb_add_si(aux2,aux2,1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_sub_si(aux2,delta,1,prec);
  arb_sqr(aux2,aux2,prec);
  arb_sub_si(aux2,aux2,4,prec);
  arb_neg(aux2,aux2);
  arb_div(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  // First line
  
  arb_set_str(aux,"2",prec);
  arb_sub(aux,aux,delta,prec);
  arb_sqrt(aux,aux,prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"1",prec);
  arb_sub(aux,aux,delta,prec);
  arb_sqrt(aux,aux,prec);
  arb_add(out,out,aux,prec);
  arb_sqrt(aux,delta,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  arb_sub_si(out,out,1,prec);
  arb_set_str(aux,"2",prec);
  arb_sqrt(aux,aux,prec);
  arb_sub(out,out,aux,prec);
  
  arb_sqrt(aux,Ips->x,prec);
  arb_mul(out,out,aux,prec);
  arb_div(out,out,Ips->u,prec);
  arb_set_si(aux,2);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(aux,"1",prec);
  arb_sub(aux,aux,delta,prec);
  arb_log(aux,aux,prec);
  arb_abs(aux,aux);
  arb_mul(out,out,aux,prec);
     
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  return 1;
}


int compute_C_Delta(arb_ptr out, const arb_t delta, slong prec){
  arb_t aux;
  arb_init(aux);
  arb_add_si(aux,delta,1,prec);
  arb_sqr(aux,aux,prec);
  arb_sub_si(aux,aux,4,prec);
  arb_neg(aux,aux);
  arb_sqr(out,delta,prec);
  arb_add_si(out,out,1,prec);
  arb_mul_si(out,out,2,prec);
  arb_div(out,out,aux,prec);
  arb_clear(aux);
  return 1;
}

int compute_c_delta(arb_ptr out, const arb_t delta, slong prec){
  arb_t aux;
  arb_init(aux);
  arb_add_si(out,delta,1,prec);
  arb_sqrt(out,out,prec);
  arb_inv(out,out,prec);
  arb_sub_si(aux,delta,1,prec);
  arb_neg(aux,aux);
  arb_sqrt(aux,aux,prec);
  arb_inv(aux,aux,prec);
  arb_add(out,out,aux,prec);
  arb_clear(aux);
  return 1;
}
