/* constant_functions.c: coefficients, remainders and scalar constants of the
 * cusp expansions of u_0 and L u_0 (Section 5 and Appendix B of the paper).
 * See constant_functions.h for the meaning of each routine.
 */
#include "constant_functions.h"

int compute_a_jk_0(arb_ptr out, int j, int k, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);
  
  if (j == 0){
    if (k < 0 || k >= ups-> N0){
      arb_clear(argument);
      arb_clear(aux);
      return 0;
    }
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);
    
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_div_si(argument,argument,2,prec);
    arb_sin_pi(aux,argument,prec);
    
    arb_mul(out,out,aux,prec);
    arb_mul(out,out,ups->a0k+k,prec);
    arb_neg(out,out);
    
    arb_clear(argument);
    arb_clear(aux);
    return 1;
  }
  else if (j == 1){
    if (k < 0 || k >= ups-> N1){
      arb_clear(argument);
      arb_clear(aux);
      return 0;
    }
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);
    
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_div_si(argument,argument,2,prec);
    arb_sin_pi(aux,argument,prec);
    
    arb_mul(out,out,aux,prec);
    arb_mul(out,out,ups->a1k+k,prec);
    arb_neg(out,out);
    
    arb_clear(argument);
    arb_clear(aux);
    return 1;
  }
  else {
    arb_clear(argument);
    arb_clear(aux);
    return 0;
  }
  
}




int compute_a_m_0(arb_ptr out, int m, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux;
  
  int k;
  
  arb_init(argument);
  arb_init(aux);
  arb_set_str(out,"0",prec);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a0k+k,prec);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups -> N1; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_add(out,out,aux,prec);
  }
  
  if (m%2 == 0){
    arb_neg(out,out);
  }
  arb_fac_ui(aux,2*m,prec);
  arb_div(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 1;
  
}



int compute_A_jk_0(arb_ptr out, int j, int k, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux;
  
  arb_init(argument);
  arb_init(aux);
  
  if (j == 0){
    if (k < 0 || k >= ups-> N0){
      arb_clear(argument);
      arb_clear(aux);
      return 0;
    }
    
    if (k == 0){ // The limit is pi/2*a_{00}
      arb_const_pi(out,prec);
      arb_div_si(out,out,2,prec);
      arb_mul(out,out,ups->a0k,prec);
      arb_clear(argument);
      arb_clear(aux);
      return 1;
    }
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"1.0",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);  
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_div_si(argument,argument,2,prec);
    arb_sin_pi(aux,argument,prec);
    
    arb_mul(out,out,aux,prec);
    arb_mul(out,out,ups->a0k+k,prec);
    
    arb_clear(argument);
    arb_clear(aux);
    return 1;
  }
  else if (j == 1){
    if (k < 0 || k >= ups-> N1){
      arb_clear(argument);
      arb_clear(aux);
      return 0;
    }
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"1.0",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);
    
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_div_si(argument,argument,2,prec);
    arb_sin_pi(aux,argument,prec);
    
    arb_mul(out,out,aux,prec);
    arb_mul(out,out,ups->a1k+k,prec);
    
    arb_clear(argument);
    arb_clear(aux);
    return 1;
  }
  else {
    arb_clear(argument);
    arb_clear(aux);
    return 0;
  }
  
}



int compute_A_m_0(arb_ptr out, int m, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux;
  
  int k;
  
  arb_init(argument);
  arb_init(aux);
  arb_set_str(out,"0",prec);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"2.0",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a0k+k,prec);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups -> N1; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"2.0",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_add(out,out,aux,prec);
  }
  
  if (m%2 == 0){
    arb_neg(out,out);
  }
  arb_fac_ui(aux,2*m,prec);
  arb_div(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 1;
  
}



int compute_b_0_0(arb_ptr out, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  int n;
  
  arb_init(aux);
  compute_a_m_0(out,1,params,prec);
  
  for (n = 1; n <= ups->N2; n++){
    arb_mul_si(aux,ups->bi+n-1,n*n,prec);
    arb_div_si(aux,aux,2,prec);
    arb_sub(out,out,aux,prec);
  }

  arb_clear(aux);
  return 1;
  
}


int compute_b_1_0(arb_ptr out, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  int n;
  
  arb_init(aux);
  compute_a_m_0(out,2,params,prec);
  
  for (n = 1; n <= ups->N2; n++){
    arb_mul_si(aux,ups->bi+n-1,n*n*n*n,prec);
    arb_div_si(aux,aux,24,prec);
    arb_add(out,out,aux,prec);
  }

  arb_clear(aux);
  return 1;
  
}

int compute_E_x_u0(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux, aux2;
  arb_t twopi;
  
  int k, m;
  
  arb_init(argument);
  arb_init(aux);
  arb_init(aux2);
  arb_init(twopi);
  arb_set_str(out,"0",prec);
  arb_const_pi(twopi,prec);
  arb_mul_si(twopi,twopi,2,prec);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_pow(aux2,twopi,argument,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_sub(argument,aux,argument,prec);
    arb_add_si(argument,argument,2*M,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a0k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_pow(aux2,twopi,argument,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_sub(argument,aux,argument,prec);
    arb_add_si(argument,argument,2*M,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  // we multiply by 2*x^(2*M)/(4*pi^2 - x^2)*(2*pi)^(5/2-2*M)
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  arb_pow_ui(aux,inp,2*M,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(argument,"2.5",prec);
  arb_sub_ui(argument,argument,2*M,prec);
  arb_pow(aux,twopi,argument,prec);
  arb_mul(out,out,aux,prec);
  arb_mul_si(out,out,2,prec);
  
  // second line
  for (m = 1; m <= ups->N2; m++){
    arb_mul_si(aux,inp,m,prec);
    arb_pow_ui(aux,aux,6,prec);
    arb_mul(aux,aux,ups->bi+m-1,prec);
    arb_abs(aux,aux);
    arb_div_si(aux,aux,720,prec);
    arb_add(out,out,aux,prec);
  }
  
  for (m = 3; m <= M-1; m++){
    compute_a_m_0(aux,m,params,prec);
    arb_abs(aux,aux);
    arb_pow_ui(aux2,inp,2*m,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  arb_clear(argument);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(twopi);
  return 1;
  
}


int compute_E_x_Lu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux, aux2;
  arb_t twopi;
  
  int k, m;
  
  arb_init(argument);
  arb_init(aux);
  arb_init(aux2);
  arb_init(twopi);
  arb_set_str(out,"0",prec);
  arb_const_pi(twopi,prec);
  arb_mul_si(twopi,twopi,2,prec);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_pow(aux2,twopi,argument,prec);
    arb_neg(argument,argument);
    arb_add_si(argument,argument,2*M-1,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a0k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_pow(aux2,twopi,argument,prec);
    arb_neg(argument,argument);
    arb_add_si(argument,argument,2*M-1,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  // we multiply by 2*x^(2*M)/(4*pi^2 - x^2)*(2*pi)^(3-2*M)
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  arb_pow_ui(aux,inp,2*M,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(argument,"3",prec);
  arb_sub_si(argument,argument,2*M,prec);
  arb_pow(aux,twopi,argument,prec);
  arb_mul(out,out,aux,prec);
  arb_mul_si(out,out,2,prec);
  
  // second line
  for (m = 1; m <= ups->N2; m++){
    arb_set_si(aux,m);
    arb_set_str(aux2,"1.5",prec);
    arb_pow(aux,aux,aux2,prec);
    
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_mul(aux,aux,ups->bi+m-1,prec);
    arb_abs(aux,aux);
    arb_div_si(aux,aux,720,prec);
    
    arb_pow_ui(aux2,inp,6,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  for (m = 3; m <= M-1; m++){
    compute_A_m_0(aux,m,params,prec);
    arb_abs(aux,aux);
    arb_pow_ui(aux2,inp,2*m,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  arb_clear(argument);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(twopi);
  return 1;
  
}


int compute_c_hat_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  
  int k;
  
  arb_init(aux);
  arb_init(aux2);
  
  arb_zero(out);
  
  for (k = 2; k < ups->N0; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_0(aux, 0, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_add(aux,aux,ups->p1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_0(aux, 1, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  
  arb_set_str(aux,"1.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_0_0(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);


  arb_set_str(aux,"3.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_1_0(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);

  arb_set_str(aux,"-0.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_E_x_u0(aux, eps, M, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  compute_a_jk_0(aux, 0, 1, params, prec);
  arb_sub(out,out,aux,prec);

  arb_clear(aux);
  arb_clear(aux2);
  return 1;
}


int compute_c_hat_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  
  int k;
  
  arb_init(aux);
  arb_init(aux2);
  
  arb_zero(out);
  
  for (k = 2; k < ups->N0; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_0(aux, 0, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_add(aux,aux,ups->p1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_0(aux, 1, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  
  arb_set_str(aux,"1.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_0_0(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);


  arb_set_str(aux,"3.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_1_0(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);

  arb_set_str(aux,"-0.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_E_x_u0(aux, eps, M, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  arb_neg(out,out);
  
  compute_a_jk_0(aux, 0, 1, params, prec);
  arb_sub(out,out,aux,prec);

  arb_clear(aux);
  arb_clear(aux2);
  return 1;
}


int compute_c_u0_hat_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t c_hat_plus;
  arb_t aux;
  arb_t lambda;
  
  
  arb_init(c_hat_plus);
  arb_init(aux);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
    
  compute_c_hat_plus(c_hat_plus, eps, M, params, prec);
  
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_hat_plus,prec);
  arb_sub(aux,lambda,aux,prec);
  
  arb_div(out,c_hat_plus,aux,prec);

  arb_clear(c_hat_plus);
  arb_clear(aux);
  arb_clear(lambda);
  return 1;  
}


int compute_c_u0_hat_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  
  arb_t c_hat_minus;
  arb_t lambda;
  
  
  arb_init(c_hat_minus);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
    
  compute_c_hat_minus(c_hat_minus, eps, M, params, prec);
  
  arb_div(out,c_hat_minus,lambda,prec);

  arb_clear(c_hat_minus);
  arb_clear(lambda);
  return 1;  
}


int compute_c_p0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  compute_c_hat_plus(out,eps,M,params,prec);
  return 1;
}


int compute_cT0_1(arb_ptr out, const arb_t eps, void* params, slong prec){
  
  arb_t argument;
  arb_t aux, aux2, aux3;
  arb_t pi;
  
  int m;
  
  arb_init(argument);
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  arb_set_str(argument,"-1.5",prec);
  arb_zeta(out,argument,prec);
  arb_div(out,out,pi,prec);
  arb_div_si(out,out,2,prec);
  
  
  arb_mul(aux,pi,pi,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_mul(aux2,eps,eps,prec);
  arb_mul_ui(aux2,aux2,4,prec);
  arb_sub(aux2,aux,aux2,prec);
  arb_set_str(argument,"2.5",prec);
  arb_zeta(aux,argument,prec);
  arb_div(aux,aux,aux2,prec);
  arb_set_str(aux2,"2",prec);
  arb_sqrt(aux2,aux2,prec);
  arb_sub_si(aux2,aux2,1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,8,prec);
  arb_set_str(argument,"1.5",prec);
  arb_pow(aux2,pi,argument,prec);
  arb_div(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  
  arb_zero(aux3);
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m);
    arb_pow(aux2,aux2,argument,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_sub(aux3,aux3,aux,prec);
  }
  compute_tail_error(aux2,argument,TANH_MAX+1,prec);
  arb_add(aux,aux2,aux3,prec);
    
  arb_mul_si(aux,aux,3,prec);
  arb_div_si(aux,aux,4,prec);
  arb_div(aux,aux,pi,prec);
    
  arb_add(out,out,aux,prec);
    
  
  arb_clear(argument);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  return 1;
}


int compute_cT0_2(arb_ptr out, const arb_t eps, void* params, slong prec){
  
  arb_t argument;
  arb_t aux, aux2, aux3;
  arb_t pi;
  
  int m;
  
  arb_init(argument);
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  
  arb_mul(aux,pi,pi,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_mul(aux2,eps,eps,prec);
  arb_sub(aux2,aux,aux2,prec);
  arb_set_str(argument,"2.5",prec);
  arb_zeta(aux,argument,prec);
  arb_div(aux,aux,aux2,prec);
  arb_sqrt(aux2,pi,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux2,pi,pi,prec);
  arb_mul(aux3,eps,eps,prec);
  arb_sub(aux2,aux2,aux3,prec);
  arb_div(out,aux,aux2,prec);
  
  arb_zero(aux3);
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m);
    arb_set_str(argument,"1.5",prec);
    arb_pow(aux2,aux2,argument,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_sub(aux3,aux3,aux,prec);
  }
  compute_tail_error(aux2,argument,TANH_MAX+1,prec);
  arb_add(aux,aux2,aux3,prec);

  arb_div_si(aux,aux,4,prec);
  arb_div(aux,aux,pi,prec);
    
  arb_add(out,out,aux,prec);
    
  
  arb_clear(argument);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  return 1;
}


/* Left-hand side minus right-hand side of the local T_0 criterion (ineqT0),
   with C_B replaced by its upper bound 1 and c_{p_0} by the computed value
   c_{p_0}(eps, M); a negative result implies the inequality of the paper. */
int compute_validation_T0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;

  arb_t c_p0;
  arb_t cT01, cT02;
  arb_t aux, aux2;
  arb_t pi;
  
  arb_zero(out);
  
  arb_init(c_p0);
  arb_init(cT01); arb_init(cT02);
  arb_init(aux); arb_init(aux2);
  arb_init(pi);
  
  arb_const_pi(pi,prec);
  
  compute_c_p0(c_p0, eps, M, params, prec);
  compute_cT0_1(cT01, eps, params, prec);
  compute_cT0_2(cT02, eps, params, prec);
  
  
  arb_set_str(aux2,"0.5",prec);
  arb_sub(aux,ups->p0,aux2,prec);
  arb_pow(aux,eps,aux,prec);
  arb_mul(out,c_p0,aux,prec);
  
  arb_mul(aux,eps,eps,prec);
  arb_add(aux2,cT01,cT02,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  
  arb_set_si(aux,2);
  arb_sqrt(aux,aux,prec);
  arb_sub_si(aux,aux,2,prec);
  arb_div_si(aux,aux,2,prec);
  arb_div(aux,aux,pi,prec);
  arb_set_str(aux2,"0.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_sub(out,out,aux,prec);
  
  
  arb_clear(c_p0);
  arb_clear(cT01); arb_clear(cT02);
  arb_clear(aux); arb_clear(aux2);
  arb_clear(pi);
  return 1;
}




/* nu = 4 L u_0(0) = 8 int_0^pi K(y) u_0(y) dy, evaluated through the Fourier
   series of K and the definition of u_0: the Clausen sums give the zeta
   values, the weights sqrt(tanh n) - 1 are summed up to n = TANH_MAX and the
   tails are bounded from n = TANH_MAX + 1. The symmetric error 4 pi eps0 of
   the wave-speed proposition is added. */
int compute_nu(arb_ptr out, const arb_t eps0, void* params, slong prec){
  u_params_struct* up = params;
  
  arb_t a0kterm,a1kterm,biterm;
  arb_t aux, aux2;
  arb_t argument;
  arb_t error_tail;
  arb_t pm_one;
  
  int i,m;
  
  arb_init(a0kterm); arb_init(a1kterm); arb_init(biterm);
  arb_zero(a0kterm); arb_zero(a1kterm); arb_zero(biterm);
  
  arb_init(aux); arb_init(aux2);
  arb_init(argument);
  arb_init(error_tail); arb_init(pm_one);
  arb_zero_pm_one(pm_one);
  
  // 4 L u_0(0): main terms and finite corrections
  
  // Main term
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_zeta(aux,argument,prec);
    
    arb_set_str(aux2,"0.5",prec);
    arb_add(argument,argument,aux2,prec);
    arb_zeta(aux2,argument,prec);
    
    arb_sub(aux,aux,aux2,prec);
    arb_mul(aux,aux,up->a0k+i,prec);
    arb_add(a0kterm,a0kterm,aux,prec);
  }
  
  
  // Main term
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_zeta(aux,argument,prec);
    
    arb_set_str(aux2,"0.5",prec);
    arb_add(argument,argument,aux2,prec);
    arb_zeta(aux2,argument,prec);
    
    arb_sub(aux,aux,aux2,prec);
    arb_mul(aux,aux,up->a1k+i,prec);
    arb_add(a1kterm,a1kterm,aux,prec);
  }
  
  
  // Correction N0
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"2",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    
    compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
    arb_mul(error_tail,error_tail,up->a0k+i,prec);
    
    for (m = 1; m <= TANH_MAX; m++){
      arb_set_si(aux,m);
      arb_tanh(aux,aux,prec);
      arb_sqrt(aux,aux,prec);
      arb_sub_si(aux,aux,1,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a0k+i,prec);
      arb_sub(a0kterm,a0kterm,aux,prec);
    }
    // Tail N0
    arb_mul(error_tail,error_tail,pm_one,prec);
    arb_add(a0kterm,a0kterm,error_tail,prec);
  }

  
  // Correction N1
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"2",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    
    compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
    arb_mul(error_tail,error_tail,up->a1k+i,prec);
    
    for (m = 1; m <= TANH_MAX; m++){
      arb_set_si(aux,m);
      arb_tanh(aux,aux,prec);
      arb_sqrt(aux,aux,prec);
      arb_sub_si(aux,aux,1,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a1k+i,prec);
      arb_sub(a1kterm,a1kterm,aux,prec);
    }
    // Tail N1
    arb_mul(error_tail,error_tail,pm_one,prec);
    arb_add(a1kterm,a1kterm,error_tail,prec);
  }
  
  
  for (i = 1; i <= up->N2; i++){
    //sqrt(tanh(k)/k)-1
    arb_set_si(aux,i);
    arb_tanh(aux,aux,prec);
    arb_div_si(aux,aux,i,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    arb_mul(aux,aux,up->bi+i-1,prec);
    arb_add(biterm,biterm,aux,prec);
  }
  
  
  arb_add(out,a0kterm,a1kterm,prec);
  arb_add(out,out,biterm,prec);
  
  arb_mul_si(out,out,4,prec);
  
  // We add 4*pi*eps_0
  arb_set(aux,eps0);
  arb_const_pi(aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_mul(aux,aux,pm_one,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);
  
  arb_clear(aux); arb_clear(aux2);
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pm_one);
  return 1;
}
