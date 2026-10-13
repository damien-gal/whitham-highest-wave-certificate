/* Blah blah blah. Javi */

#include "constant_functions.h"

// Overall: 1 if OKAY, 0 if NOT OKAY

// TODO(Javi): Check if 0 is OKAY or not
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




// TODO(Javi): Check if 0 is OKAY or not
int compute_a_jk_1(arb_ptr out, int j, int k, void* params, slong prec){
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
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);
    
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_div_si(argument,argument,2,prec);
    arb_cos_pi(aux,argument,prec);
    
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
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);
    
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_div_si(argument,argument,2,prec);
    arb_cos_pi(aux,argument,prec);
    
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


// TODO(Javi): Check if 0 is OKAY or not
int compute_a_jk_2(arb_ptr out, int j, int k, void* params, slong prec){
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
    arb_set_str(aux,"-1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);
    
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
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
    arb_set_str(aux,"-1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    arb_gamma(out,argument,prec);
    
    
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
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


// TODO(Javi): Check if 0 is OKAY or not
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



// TODO(Javi): Check if 0 is OKAY or not
int compute_a_m_1(arb_ptr out, int m, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux;
  
  int k;
  
  arb_init(argument);
  arb_init(aux);
  arb_set_str(out,"0",prec);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a0k+k,prec);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups -> N1; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_add(out,out,aux,prec);
  }
  
  if (m%2 != 0){
    arb_neg(out,out);
  }
  arb_fac_ui(aux,2*m+1,prec);
  arb_div(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 1;
  
}



// TODO(Javi): Check if 0 is OKAY or not
int compute_a_m_2(arb_ptr out, int m, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux;
  
  int k;
  
  arb_init(argument);
  arb_init(aux);
  arb_set_str(out,"0",prec);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a0k+k,prec);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups -> N1; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_add(out,out,aux,prec);
  }
  
  if (m%2 != 0){
    arb_neg(out,out);
  }
  arb_fac_ui(aux,2*m,prec);
  arb_div(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 1;
  
}



// TODO(Javi): Check if 0 is OKAY or not
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



// TODO(Javi): Check if 0 is OKAY or not
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



// TODO(Javi): Check if 0 is OKAY or not
int compute_A_m_1(arb_ptr out, int m, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t argument;
  arb_t aux;
  
  int k;
  
  arb_init(argument);
  arb_init(aux);
  arb_set_str(out,"0",prec);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a0k+k,prec);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups -> N1; k++){
    arb_mul_si(argument,ups->p0,k,prec);
    arb_add(argument,argument,ups->p1,prec);
    arb_sub_si(argument,argument,2*m,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_add(out,out,aux,prec);
  }
  
  if (m%2 != 0){
    arb_neg(out,out);
  }
  arb_fac_ui(aux,2*m+1,prec);
  arb_div(out,out,aux,prec);
  
  arb_clear(argument);
  arb_clear(aux);
  return 1;
  
}


// TODO(Javi): Check if 0 is OKAY or not
int compute_A_m_2(arb_ptr out, int m, void* params, slong prec){
  compute_A_m_1(out,m,params,prec);
  arb_mul_si(out,out,2*m+1,prec);
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


int compute_b_0_1(arb_ptr out, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  int n;
  
  arb_init(aux);
  compute_a_m_1(out,0,params,prec);
  
  for (n = 1; n <= ups->N2; n++){
    arb_mul_si(aux,ups->bi+n-1,n*n,prec);
    arb_sub(out,out,aux,prec);
  }

  arb_clear(aux);
  return 1;
  
}

int compute_b_0_2(arb_ptr out, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  int n;
  
  arb_init(aux);
  compute_a_m_2(out,0,params,prec);
  
  for (n = 1; n <= ups->N2; n++){
    arb_mul_si(aux,ups->bi+n-1,n*n,prec);
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

int compute_b_1_1(arb_ptr out, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  int n;
  
  arb_init(aux);
  compute_a_m_1(out,1,params,prec);
  
  for (n = 1; n <= ups->N2; n++){
    arb_mul_si(aux,ups->bi+n-1,n*n*n*n,prec);
    arb_div_si(aux,aux,6,prec);
    arb_add(out,out,aux,prec);
  }

  arb_clear(aux);
  return 1;
  
}

int compute_b_1_2(arb_ptr out, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  int n;
  
  arb_init(aux);
  compute_a_m_2(out,1,params,prec);
  
  for (n = 1; n <= ups->N2; n++){
    arb_mul_si(aux,ups->bi+n-1,n*n*n*n,prec);
    arb_div_si(aux,aux,2,prec);
    arb_add(out,out,aux,prec);
  }

  arb_clear(aux);
  return 1;
  
}


int compute_E_x_Cz(arb_ptr out, const arb_t inp, const arb_t z, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2, auxout;
  arb_t twopi;
  
  long int M;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(auxout);
  arb_init(twopi);
  arb_set_str(auxout,"0",prec);
  arb_const_pi(twopi,prec);
  arb_mul_si(twopi,twopi,2,prec);
  
  // Calculate M
  arb_add_si(aux,z,1,prec);
  arb_div_si(aux,aux,2,prec);
  arb_ceil(aux,aux,prec);
  arb_get_unique_fmpz(&M, aux);
    
  // Calculate the bound
  arb_set_si(aux,2*M+1);
  arb_sub(aux,aux,z,prec);
  arb_zeta(auxout,aux,prec);
  
  arb_set_si(aux,-2*M+1);
  arb_add(aux,aux,z,prec);
  arb_pow(aux,twopi,aux,prec);
  arb_mul(auxout,auxout,aux,prec);
  
  arb_mul_si(auxout,auxout,2,prec);
  
  arb_pow_ui(aux,inp,2*M,prec);
  arb_mul(auxout,auxout,aux,prec);
  
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,auxout,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(auxout);
  arb_clear(twopi);
  return 1;
}


int compute_E_x_Sz(arb_ptr out, const arb_t inp, const arb_t z, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2, auxout;
  arb_t twopi;
  
  long int M;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(auxout);
  arb_init(twopi);
  arb_set_str(auxout,"0",prec);
  arb_const_pi(twopi,prec);
  arb_mul_si(twopi,twopi,2,prec);
  
  // Calculate M
  arb_add_si(aux,z,1,prec);
  arb_div_si(aux,aux,2,prec);
  arb_ceil(aux,aux,prec);
  arb_get_unique_fmpz(&M, aux);
  
  // Calculate the bound
  arb_set_si(aux,2*M+2);
  arb_sub(aux,aux,z,prec);
  arb_zeta(auxout,aux,prec);
  
  arb_set_si(aux,-2*M);
  arb_add(aux,aux,z,prec);
  arb_pow(aux,twopi,aux,prec);
  arb_mul(auxout,auxout,aux,prec);
  
  arb_mul_si(auxout,auxout,2,prec);
  
  arb_pow_ui(aux,inp,2*M+1,prec);
  arb_mul(auxout,auxout,aux,prec);
  
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,auxout,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(auxout);
  arb_clear(twopi);
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


int compute_E_x_du0(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
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
    arb_set_str(aux,"1.5",prec);
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
    arb_set_str(aux,"1.5",prec);
    arb_sub(argument,aux,argument,prec);
    arb_add_si(argument,argument,2*M,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  // we multiply by 2*x^(2*M+1)/(4*pi^2 - x^2)*(2*pi)^(1/2-2*M)
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  arb_pow_ui(aux,inp,2*M+1,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(argument,"0.5",prec);
  arb_sub_ui(argument,argument,2*M,prec);
  arb_pow(aux,twopi,argument,prec);
  arb_mul(out,out,aux,prec);
  arb_mul_si(out,out,2,prec);
  
  // second line
  for (m = 1; m <= ups->N2; m++){
    arb_mul_si(aux,inp,m,prec);
    arb_pow_ui(aux,aux,5,prec);
    arb_mul_si(aux,aux,m,prec);
    arb_mul(aux,aux,ups->bi+m-1,prec);
    arb_abs(aux,aux);
    arb_div_si(aux,aux,120,prec);
    arb_add(out,out,aux,prec);
  }
  
  for (m = 2; m <= M-1; m++){
    compute_a_m_1(aux,m,params,prec);
    arb_abs(aux,aux);
    arb_pow_ui(aux2,inp,2*m+1,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  arb_clear(argument);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(twopi);
  return 1;
  
}


int compute_E_x_ddu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
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
    arb_set_str(aux,"1.5",prec);
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
    arb_set_str(aux,"1.5",prec);
    arb_sub(argument,aux,argument,prec);
    arb_add_si(argument,argument,2*M,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  // we multiply by 2*x^(2*M)/(4*pi^2 - x^2)*(2*pi)^(1/2-2*M)
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  arb_pow_ui(aux,inp,2*M,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(argument,"0.5",prec);
  arb_sub_si(argument,argument,2*M,prec);
  arb_pow(aux,twopi,argument,prec);
  arb_mul(out,out,aux,prec);
  arb_mul_si(out,out,2,prec);
  
  // second line
  for (m = 1; m <= ups->N2; m++){
    arb_mul_si(aux,inp,m,prec);
    arb_pow_ui(aux,aux,4,prec);
    arb_mul_si(aux,aux,m,prec);
    arb_mul_si(aux,aux,m,prec);
    arb_mul(aux,aux,ups->bi+m-1,prec);
    arb_abs(aux,aux);
    arb_div_si(aux,aux,24,prec);
    arb_add(out,out,aux,prec);
  }
  
  for (m = 2; m <= M-1; m++){
    compute_a_m_2(aux,m,params,prec);
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


int compute_E_x_dLu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
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
    arb_add_si(argument,argument,2*M+1,prec);
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
    arb_add_si(argument,argument,2*M+1,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  // we multiply by 2*x^(2*M+1)/(4*pi^2 - x^2)*(2*pi)^(1-2*M)
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  arb_pow_ui(aux,inp,2*M+1,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(argument,"1",prec);
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
    arb_div_si(aux,aux,120,prec);
    
    arb_pow_ui(aux2,inp,5,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  for (m = 2; m <= M-1; m++){
    compute_A_m_1(aux,m,params,prec);
    arb_abs(aux,aux);
    arb_pow_ui(aux2,inp,2*m+1,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  arb_clear(argument);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(twopi);
  return 1;
  
}


int compute_E_x_ddLu0(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
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
    arb_add_si(argument,argument,2*M+1,prec);
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
    arb_add_si(argument,argument,2*M+1,prec);
    arb_zeta(aux,argument,prec);
    arb_mul(aux,aux,ups->a1k+k,prec);
    arb_abs(aux,aux);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  // we multiply by 2*x^(2*M)/(4*pi^2 - x^2)*(2*pi)^(1-2*M)
  arb_sqr(aux,twopi,prec);
  arb_sqr(aux2,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  arb_pow_ui(aux,inp,2*M,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(argument,"1",prec);
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
    arb_div_si(aux,aux,24,prec);
    
    arb_pow_ui(aux2,inp,4,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  for (m = 2; m <= M-1; m++){
    compute_A_m_2(aux,m,params,prec);
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


int compute_E_x_f(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t c_ep;
  arb_t pi;
  arb_t lambda;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(c_ep);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux,aux2,1,prec);
    
    arb_set_si(aux2,m);
    arb_sqrt(aux2,aux2,prec);
    
    arb_div(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  arb_set_str(aux2,"-0.5",prec);
  compute_tail_error(aux,aux2,TANH_MAX+1,prec);
  arb_add(out,out,aux,prec);
  
  compute_cp_p0_plus(c_ep,inp,M,params,prec);
  arb_mul_si(aux,c_ep,2,prec);
  compute_cpp_p0_plus(c_ep,inp,M,params,prec);  
  arb_add(aux,aux,c_ep,prec);
  arb_div(aux,aux,lambda,prec);
  arb_pow(aux2,inp,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_str(aux2,"19",prec);
  arb_div_si(aux2,aux2,8,prec);
  arb_add(aux,aux,aux2,prec);
  arb_mul(out,out,aux,prec);
  
  arb_mul_si(out,out,2,prec);
  arb_div(out,out,pi,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(c_ep);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}


int compute_E_x_F0(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2, aux3;
  
  int k,j,n;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_set_str(out,"0",prec);
  
  // First line
  compute_A_jk_0(out,1,0,params,prec);
  compute_a_jk_0(aux,0,0,params,prec);
  compute_a_jk_0(aux2,1,0,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  arb_abs(out,out);
  arb_sub(aux,ups->p1,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  arb_mul(out,out,aux,prec);
  
  
  compute_A_jk_0(aux,1,1,params,prec);
  compute_a_jk_0(aux2,0,0,params,prec);
  compute_a_jk_0(aux3,1,1,params,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_sub(aux,aux,aux2,prec);
  compute_a_jk_0(aux2,0,1,params,prec);
  compute_a_jk_0(aux3,1,0,params,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_pow(aux2,inp,ups->p1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  // Second line
  for (k = 2; k < ups->N0; k++){
    compute_A_jk_0(aux,0,k,params,prec);
    if (k%2 == 0){
      compute_a_jk_0(aux2,0,k/2,params,prec);
      arb_sqr(aux2,aux2,prec);
      arb_sub(aux,aux,aux2,prec);
    }
    for (j = 0; j <= (k-1)/2; j++){
      compute_a_jk_0(aux2,0,j,params,prec);
      compute_a_jk_0(aux3,0,k-j,params,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul_si(aux2,aux2,2,prec);
      arb_sub(aux,aux,aux2,prec);
    }
    arb_abs(aux,aux);
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_pow(aux2,inp,aux2,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  // Third line
  
  compute_A_m_0(aux, 1, params, prec);

  for (k = 1; k <= ups->N2; k++){
    arb_set_si(aux2,k);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,k*k*k,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_mul(aux2,aux2,ups->bi+k-1,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_sub(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N0; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a0k+k,prec);
      arb_div_si(aux2,aux2,2,prec);
      arb_add(aux,aux,aux2,prec);
    }
  }
  
  for (k = 0; k < ups->N1; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_sub(aux2,aux2,ups->p1,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a1k+k,prec);
      arb_div_si(aux2,aux2,2,prec);
      arb_add(aux,aux,aux2,prec);
    }
  }
  
  arb_abs(aux,aux);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a0k+k,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_sub(aux2,aux2,ups->p1,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a1k+k,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  arb_set_si(aux2,1);
  arb_sub(aux2,aux2,ups->p0,prec);
  arb_pow(aux2,inp,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
      
  // Fourth line
  
  compute_A_m_0(aux, 2, params, prec);
  
  for (k = 1; k <= ups->N2; k++){
    arb_set_si(aux2,k);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,k*k*k*k*k*k*k,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_mul(aux2,aux2,ups->bi+k-1,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N0; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_add_si(aux2,aux2,2,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a0k+k,prec);
      arb_div_si(aux2,aux2,24,prec);
      arb_sub(aux,aux,aux2,prec);
    }
  }
  
  for (k = 0; k < ups->N1; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_sub(aux2,aux2,ups->p1,prec);
      arb_add_si(aux2,aux2,2,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a1k+k,prec);
      arb_div_si(aux2,aux2,24,prec);
      arb_sub(aux,aux,aux2,prec);
    }
  }
  
  compute_b_0_0(aux2,params,prec);
  arb_mul(aux2,aux2,aux2,prec);
  arb_sub(aux,aux,aux2,prec);
  
  arb_abs(aux,aux);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_add_si(aux2,aux2,2,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a0k+k,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_sub(aux2,aux2,ups->p1,prec);
    arb_add_si(aux2,aux2,2,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a1k+k,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  arb_set_si(aux2,3);
  arb_sub(aux2,aux2,ups->p0,prec);
  arb_pow(aux2,inp,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  // Fifth line
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"1.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"1.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_b_0_0(aux2,params,prec);
  
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"3.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"3.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_b_1_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  arb_set_si(aux,7);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  compute_b_1_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  // Sixth line
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_E_x_u0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);

  arb_set_si(aux,-1);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  compute_E_x_Lu0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);

  arb_set_si(aux,-1);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  compute_E_x_u0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  return 1;
}


int compute_E_x_F1(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2, aux3;
  
  int k,j,n;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_set_str(out,"0",prec);
  
  // First line
  compute_A_jk_0(out,1,0,params,prec);
  compute_a_jk_0(aux,0,0,params,prec);
  compute_a_jk_0(aux2,1,0,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  arb_abs(out,out);
  arb_sub(aux,ups->p1,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_add_si(aux,ups->p1,1,prec);
  arb_mul(out,out,aux,prec);
  
  
  
  compute_A_jk_0(aux,1,1,params,prec);
  compute_a_jk_0(aux2,0,0,params,prec);
  compute_a_jk_0(aux3,1,1,params,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_sub(aux,aux,aux2,prec);
  compute_a_jk_0(aux2,0,1,params,prec);
  compute_a_jk_0(aux3,1,0,params,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_pow(aux2,inp,ups->p1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(aux2,ups->p0,ups->p1,prec);
  arb_add_si(aux2,aux2,1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  // Second line
  for (k = 2; k < ups->N0; k++){
    compute_A_jk_0(aux,0,k,params,prec);
    if (k%2 == 0){
      compute_a_jk_0(aux2,0,k/2,params,prec);
      arb_sqr(aux2,aux2,prec);
      arb_sub(aux,aux,aux2,prec);
    }
    for (j = 0; j <= (k-1)/2; j++){
      compute_a_jk_0(aux2,0,j,params,prec);
      compute_a_jk_0(aux3,0,k-j,params,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul_si(aux2,aux2,2,prec);
      arb_sub(aux,aux,aux2,prec);
    }
    arb_abs(aux,aux);
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_pow(aux2,inp,aux2,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_mul_si(aux2,ups->p0,k,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  // Third line
  
  compute_A_m_0(aux, 1, params, prec);
  
  for (k = 1; k <= ups->N2; k++){
    arb_set_si(aux2,k);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,k*k*k,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_mul(aux2,aux2,ups->bi+k-1,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_sub(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N0; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a0k+k,prec);
      arb_div_si(aux2,aux2,2,prec);
      arb_add(aux,aux,aux2,prec);
    }
  }
  
  for (k = 0; k < ups->N1; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_sub(aux2,aux2,ups->p1,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a1k+k,prec);
      arb_div_si(aux2,aux2,2,prec);
      arb_add(aux,aux,aux2,prec);
    }
  }
  
  arb_abs(aux,aux);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a0k+k,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_sub(aux2,aux2,ups->p1,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a1k+k,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  arb_set_si(aux2,1);
  arb_sub(aux2,aux2,ups->p0,prec);
  arb_pow(aux2,inp,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
      
  // Fourth line
  
  compute_A_m_0(aux, 2, params, prec);
  
  for (k = 1; k <= ups->N2; k++){
    arb_set_si(aux2,k);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,k*k*k*k*k*k*k,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_mul(aux2,aux2,ups->bi+k-1,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N0; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_add_si(aux2,aux2,2,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a0k+k,prec);
      arb_div_si(aux2,aux2,24,prec);
      arb_sub(aux,aux,aux2,prec);
    }
  }
  
  for (k = 0; k < ups->N1; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_sub(aux2,aux2,ups->p1,prec);
      arb_add_si(aux2,aux2,2,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a1k+k,prec);
      arb_div_si(aux2,aux2,24,prec);
      arb_sub(aux,aux,aux2,prec);
    }
  }
  
  compute_b_0_0(aux2,params,prec);
  arb_mul(aux2,aux2,aux2,prec);
  arb_sub(aux,aux,aux2,prec);
  
  arb_abs(aux,aux);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_add_si(aux2,aux2,2,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a0k+k,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_sub(aux2,aux2,ups->p1,prec);
    arb_add_si(aux2,aux2,2,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a1k+k,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  arb_set_si(aux2,3);
  arb_sub(aux2,aux2,ups->p0,prec);
  arb_pow(aux2,inp,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_add(out,out,aux,prec);
  
  // Fifth line
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"1.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"1.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_b_0_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"3.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"3.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_b_1_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  arb_set_si(aux,7);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  compute_b_1_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,8,prec);
  arb_add(out,out,aux,prec);
  
  // Sixth line
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_E_x_u0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_E_x_du0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);

  arb_neg(aux,ups->p0);
  arb_pow(aux,inp,aux,prec);
  compute_E_x_dLu0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);

  arb_neg(aux,ups->p0);
  arb_pow(aux,inp,aux,prec);
  compute_E_x_u0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  compute_E_x_du0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  return 1;
}


int compute_E_x_F2(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2, aux3;
  
  int k,j,n;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_set_str(out,"0",prec);
  
  // First line
  compute_A_jk_0(out,1,0,params,prec);
  compute_a_jk_0(aux,0,0,params,prec);
  compute_a_jk_0(aux2,1,0,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  arb_abs(out,out);
  arb_sub(aux,ups->p1,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_add_si(aux,ups->p1,1,prec);
  arb_mul(out,out,aux,prec);
  arb_mul(out,out,ups->p1,prec);
  
  
  
  compute_A_jk_0(aux,1,1,params,prec);
  compute_a_jk_0(aux2,0,0,params,prec);
  compute_a_jk_0(aux3,1,1,params,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_sub(aux,aux,aux2,prec);
  compute_a_jk_0(aux2,0,1,params,prec);
  compute_a_jk_0(aux3,1,0,params,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_pow(aux2,inp,ups->p1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(aux2,ups->p0,ups->p1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add_si(aux2,aux2,1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  // Second line
  for (k = 2; k < ups->N0; k++){
    compute_A_jk_0(aux,0,k,params,prec);
    if (k%2 == 0){
      compute_a_jk_0(aux2,0,k/2,params,prec);
      arb_sqr(aux2,aux2,prec);
      arb_sub(aux,aux,aux2,prec);
    }
    for (j = 0; j <= (k-1)/2; j++){
      compute_a_jk_0(aux2,0,j,params,prec);
      compute_a_jk_0(aux3,0,k-j,params,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul_si(aux2,aux2,2,prec);
      arb_sub(aux,aux,aux2,prec);
    }
    arb_abs(aux,aux);
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_pow(aux2,inp,aux2,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_mul_si(aux2,ups->p0,k,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  // Third line
  
  compute_A_m_0(aux, 1, params, prec);
  
  for (k = 1; k <= ups->N2; k++){
    arb_set_si(aux2,k);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,k*k*k,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_mul(aux2,aux2,ups->bi+k-1,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_sub(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N0; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a0k+k,prec);
      arb_div_si(aux2,aux2,2,prec);
      arb_add(aux,aux,aux2,prec);
    }
  }
  
  for (k = 0; k < ups->N1; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_sub(aux2,aux2,ups->p1,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a1k+k,prec);
      arb_div_si(aux2,aux2,2,prec);
      arb_add(aux,aux,aux2,prec);
    }
  }
  
  arb_abs(aux,aux);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a0k+k,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_sub(aux2,aux2,ups->p1,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a1k+k,prec);
    arb_div_si(aux2,aux2,2,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  arb_set_si(aux2,1);
  arb_sub(aux2,aux2,ups->p0,prec);
  arb_pow(aux2,inp,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
      
  // Fourth line
  
  compute_A_m_0(aux, 2, params, prec);
  
  for (k = 1; k <= ups->N2; k++){
    arb_set_si(aux2,k);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,k*k*k*k*k*k*k,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_mul(aux2,aux2,ups->bi+k-1,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N0; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_add_si(aux2,aux2,2,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a0k+k,prec);
      arb_div_si(aux2,aux2,24,prec);
      arb_sub(aux,aux,aux2,prec);
    }
  }
  
  for (k = 0; k < ups->N1; k++){
    for (n = 1; n <= TANH_MAX; n++){
      arb_mul_si(aux2,ups->p0,-k,prec);
      arb_sub(aux2,aux2,ups->p1,prec);
      arb_add_si(aux2,aux2,2,prec);
      arb_set_si(aux3,n);
      arb_pow(aux2,aux3,aux2,prec);
      arb_set_si(aux3,n);
      arb_tanh(aux3,aux3,prec);
      arb_sqrt(aux3,aux3,prec);
      arb_sub_si(aux3,aux3,1,prec);
      arb_mul(aux2,aux2,aux3,prec);
      arb_mul(aux2,aux2,ups->a1k+k,prec);
      arb_div_si(aux2,aux2,24,prec);
      arb_sub(aux,aux,aux2,prec);
    }
  }
  
  compute_b_0_0(aux2,params,prec);
  arb_mul(aux2,aux2,aux2,prec);
  arb_sub(aux,aux,aux2,prec);
  
  arb_abs(aux,aux);
  
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_add_si(aux2,aux2,2,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a0k+k,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,-k,prec);
    arb_sub(aux2,aux2,ups->p1,prec);
    arb_add_si(aux2,aux2,2,prec);
    compute_tail_error(aux2,aux2,TANH_MAX+1,prec);
    arb_mul(aux2,aux2,ups->a1k+k,prec);
    arb_div_si(aux2,aux2,24,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  arb_set_si(aux2,3);
  arb_sub(aux2,aux2,ups->p0,prec);
  arb_pow(aux2,inp,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,12,prec);
  arb_add(out,out,aux,prec);
  
  // Fifth line
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"1.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_mul(aux3,aux2,aux3,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"1.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_mul(aux3,aux2,aux3,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_b_0_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"3.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_mul(aux3,aux2,aux3,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"3.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_mul(aux3,aux2,aux3,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_b_1_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  arb_set_si(aux,7);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  compute_b_1_0(aux2,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,56,prec);
  arb_add(out,out,aux,prec);
  
  // Sixth line
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_mul(aux3,aux2,aux3,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_mul(aux3,aux2,aux3,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_E_x_u0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"-0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux3,inp,aux2,prec);
    arb_add(aux2,aux2,ups->p0,prec);
    arb_add_si(aux2,aux2,1,prec);
    arb_mul(aux2,aux2,aux3,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_E_x_du0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  
  arb_set_str(aux,"0",prec);
  for (k = 0; k < ups->N0; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_set_str(aux3,"0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,0,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux2,ups->p0,k-1,prec);
    arb_add(aux2,aux2,ups->p1,prec);
    arb_set_str(aux3,"0.5",prec);
    arb_add(aux2,aux2,aux3,prec);
    arb_pow(aux2,inp,aux2,prec);
    compute_a_jk_0(aux3,1,k,params,prec);
    arb_mul(aux2,aux2,aux3,prec);
    arb_abs(aux2,aux2);
    arb_add(aux,aux,aux2,prec);
  }
  
  compute_E_x_ddu0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  

  arb_neg(aux,ups->p0);
  arb_add_si(aux,aux,1,prec);
  arb_pow(aux,inp,aux,prec);
  compute_E_x_ddLu0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);

  arb_neg(aux,ups->p0);
  arb_add_si(aux,aux,1,prec);
  arb_pow(aux,inp,aux,prec);
  compute_E_x_du0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);

  arb_neg(aux,ups->p0);
  arb_add_si(aux,aux,1,prec);
  arb_pow(aux,inp,aux,prec);
  compute_E_x_u0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  compute_E_x_ddu0(aux2,inp,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
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
  u_params_struct* ups = params;
  
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


int compute_cp_p0_zero(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  
  int k;
  
  arb_init(aux);
  arb_init(aux2);
  
  arb_zero(out);
  
  for (k = 1; k < ups->N0; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_1(aux, 0, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_add(aux,aux,ups->p1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_1(aux, 1, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  
  arb_set_str(aux,"1.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_0_1(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);


  arb_set_str(aux,"3.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_1_1(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);

  arb_set_str(aux,"0.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_E_x_du0(aux, eps, M, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  return 1;
}


int compute_cp_p0_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  arb_t c_ep;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(c_ep);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_cp_p0_zero(out,eps,M,params,prec);
  compute_c_u0_hat_plus(c_ep,eps,M,params,prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_add_si(aux,aux,1,prec);
  arb_mul(out,out,aux,prec);
  
  arb_mul(aux,lambda,c_ep,prec);
  arb_div_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec); 
  
  arb_clear(aux);
  arb_clear(c_ep);
  arb_clear(lambda);
  return 1;

}

int compute_cp_p0_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  arb_t c_ep;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(c_ep);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_cp_p0_zero(out,eps,M,params,prec);
  compute_c_u0_hat_plus(c_ep,eps,M,params,prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_add_si(aux,aux,1,prec);
  arb_mul(out,out,aux,prec);
  
  compute_c_u0_hat_minus(c_ep,eps,M,params,prec);
  arb_mul(aux,lambda,c_ep,prec);
  arb_div_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec); 
  
  arb_clear(aux);
  arb_clear(c_ep);
  arb_clear(lambda);
  return 1;

}


int compute_cpp_p0_zero(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  
  int k;
  
  arb_init(aux);
  arb_init(aux2);
  
  arb_zero(out);
  
  for (k = 1; k < ups->N0; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_2(aux, 0, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  
  for (k = 0; k < ups->N1; k++){
    arb_mul_si(aux,ups->p0,k-1,prec);
    arb_add(aux,aux,ups->p1,prec);
    arb_pow(aux2,eps,aux,prec);
    compute_a_jk_2(aux, 1, k, params, prec);
    arb_mul(aux,aux,aux2,prec);
    arb_abs(aux,aux);
    arb_add(out,out,aux,prec);
  }
  
  arb_set_str(aux,"1.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_0_2(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);


  arb_set_str(aux,"3.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_b_1_2(aux,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);

  arb_set_str(aux,"1.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux2,eps,aux,prec);
  compute_E_x_ddu0(aux, eps, M, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_abs(aux,aux);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  return 1;
}


int compute_cpp_p0_plus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  arb_t c_ep;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(c_ep);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_cpp_p0_zero(out,eps,M,params,prec);
  compute_c_u0_hat_plus(c_ep,eps,M,params,prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_add_si(aux,aux,1,prec);
  arb_mul(out,out,aux,prec);
  
  arb_mul(aux,lambda,c_ep,prec);
  arb_div_si(aux,aux,4,prec);
  arb_add(out,out,aux,prec);
  
  compute_cp_p0_plus(c_ep,eps,M,params,prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_div(aux,aux,lambda,prec);
  arb_add_si(aux,aux,1,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  
  arb_clear(aux);
  arb_clear(c_ep);
  arb_clear(lambda);
  return 1;
}

int compute_cpp_p0_minus(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  arb_t c_ep;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(c_ep);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_cpp_p0_zero(out,eps,M,params,prec);
  compute_c_u0_hat_minus(c_ep,eps,M,params,prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_add_si(aux,aux,1,prec);
  arb_mul(out,out,aux,prec);
  
  arb_mul(aux,lambda,c_ep,prec);
  arb_div_si(aux,aux,4,prec);
  arb_sub(out,out,aux,prec);
  
  compute_cp_p0_plus(c_ep,eps,M,params,prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_div(aux,aux,lambda,prec);
  arb_add_si(aux,aux,1,prec);
  compute_cp_p0_minus(c_ep,eps,M,params,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(c_ep);
  arb_clear(lambda);
  return 1;
}



int compute_c_F0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_A_jk_0(out, 0, 1, params, prec);
  
  compute_a_jk_0(aux, 0, 0, params, prec);
  compute_a_jk_0(aux2, 0, 1, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  arb_abs(out,out);
  compute_E_x_F0(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"0.5",prec);
  arb_sub(aux,ups->p0,aux,prec);
  arb_pow(aux,eps,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,lambda,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}


int compute_c_F1(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_A_jk_0(out, 0, 1, params, prec);
  
  compute_a_jk_0(aux, 0, 0, params, prec);
  compute_a_jk_0(aux2, 0, 1, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  arb_abs(out,out);
  arb_add_si(aux,ups->p0,1,prec);
  arb_mul(out,out,aux,prec);
  compute_E_x_F1(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"0.5",prec);
  arb_sub(aux,ups->p0,aux,prec);
  arb_pow(aux,eps,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,lambda,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}

int compute_c_F2(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_A_jk_0(out, 0, 1, params, prec);
  
  compute_a_jk_0(aux, 0, 0, params, prec);
  compute_a_jk_0(aux2, 0, 1, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sub(out,out,aux,prec);
  arb_abs(out,out);
  arb_add_si(aux,ups->p0,1,prec);
  arb_mul(out,out,aux,prec);
  arb_mul(out,out,ups->p0,prec);
  compute_E_x_F2(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"0.5",prec);
  arb_sub(aux,ups->p0,aux,prec);
  arb_pow(aux,eps,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,lambda,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}


int compute_bound_alpha_eps_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t c_ep;
  arb_t aux, aux2;
  arb_t lambda;
  
  arb_init(c_ep);
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
    
  compute_c_u0_hat_plus(c_ep, eps, M, params, prec);
  
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_add_si(aux,aux,1,prec);
  
  arb_sqrt(aux2,eps,prec);
  arb_mul(aux,aux,aux2,prec);
  
  arb_div(aux,aux,lambda,prec);
  arb_div_si(out,aux,2,prec);

  arb_clear(c_ep);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}


int compute_bound_alpha1(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t c_ep;
  arb_t alpha0;
  arb_t aux, aux2;
  arb_t lambda;
  
  
  arb_init(c_ep);
  arb_init(alpha0);
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  compute_bound_alpha_eps_0(alpha0, eps, M, params, prec);
    
  compute_cp_p0_plus(c_ep, eps, M, params, prec);
  
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_div(aux,aux,lambda,prec);
  
  arb_set_str(aux2,"0.5",prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_mul(out,aux,alpha0,prec);

  arb_clear(c_ep);
  arb_clear(alpha0);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}

int compute_bound_alpha2(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t c_ep;
  arb_t alpha0;
  arb_t aux, aux2;
  arb_t lambda;
  
  
  arb_init(c_ep);
  arb_init(alpha0);
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  compute_bound_alpha_eps_0(alpha0, eps, M, params, prec);
    
  compute_cpp_p0_plus(c_ep, eps, M, params, prec);
  
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(aux,aux,c_ep,prec);
  arb_div(aux,aux,lambda,prec);
  
  arb_set_str(aux2,"0.75",prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_mul(out,aux,alpha0,prec);

  arb_clear(c_ep);
  arb_clear(alpha0);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}


int compute_bound_alpha2bar(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  compute_bound_alpha1(out,eps,M,params,prec);
  arb_mul_si(out,out,4,prec);
  return 1;
}


int compute_bound_alphaf(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t c_ep;
  arb_t aux, aux2, aux3;
  arb_t lambda;
  arb_t twopi;
  
  arb_init(c_ep);
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(lambda);
  arb_init(twopi);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  arb_const_pi(twopi,prec);
  arb_mul_si(twopi,twopi,2,prec);
  
  // Third line
  compute_E_x_f(out, eps, M, params, prec);
  arb_mul(out,out,twopi,prec);
  arb_div_si(out,out,4,prec);
  
  arb_set_str(aux,"2.5",prec);
  compute_E_x_Cz(aux, eps, aux, params, prec);
  arb_mul_si(aux,aux,3,prec);
  arb_div_si(aux,aux,4,prec);
  arb_div(aux,aux,eps,prec);
  arb_div(aux,aux,eps,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"1.5",prec);
  compute_E_x_Sz(aux, eps, aux, params, prec);
  arb_div(aux,aux,eps,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"0.5",prec);
  compute_E_x_Cz(aux, eps, aux, params, prec);
  arb_add(out,out,aux,prec);
  
  // Second line
  
  arb_set_str(aux2,"2.5",prec);
  compute_E_x_Cz(aux2, eps, aux2, params, prec);
  arb_div(aux2,aux2,eps,prec);
  arb_div(aux2,aux2,eps,prec);
  
  arb_set_str(aux3,"0.5",prec);
  arb_zeta(aux3,aux3,prec);
  arb_div_si(aux3,aux3,2,prec);
  arb_add(aux2,aux2,aux3,prec);
  compute_cpp_p0_plus(c_ep,eps,M,params,prec);
  arb_mul(aux,aux2,c_ep,prec);
  
  
  arb_set_str(aux2,"1.5",prec);
  compute_E_x_Sz(aux2, eps, aux2, params, prec);
  arb_div(aux2,aux2,eps,prec);
  
  arb_set_str(aux3,"-1.5",prec);
  arb_zeta(aux3,aux3,prec);
  arb_div_si(aux3,aux3,6,prec);
  arb_mul(aux3,aux3,eps,prec);
  arb_mul(aux3,aux3,eps,prec);
  arb_sub(aux2,aux2,aux3,prec);
  
  arb_set_str(aux3,"0.5",prec);
  arb_zeta(aux3,aux3,prec);
  arb_add(aux2,aux2,aux3,prec);
  
  compute_cp_p0_plus(c_ep,eps,M,params,prec);
  arb_mul(aux2,aux2,c_ep,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_add(aux,aux,aux2,prec);
  
  arb_div(aux,aux,lambda,prec);
  arb_pow(aux2,eps,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  
  arb_add(out,out,aux,prec);
  
  // First line
  compute_cpp_p0_plus(c_ep, eps, M, params, prec);
  arb_div_si(aux,c_ep,3,prec);
  compute_cp_p0_plus(c_ep, eps, M, params, prec);
  arb_add(aux,aux,c_ep,prec);
  arb_sqrt(aux2,twopi,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_div(aux,aux,lambda,prec);
  arb_set_str(aux2,"0.5",prec);
  arb_sub(aux2,ups->p0,aux2,prec);
  arb_pow(aux2,eps,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  
  arb_set_str(aux2,"0.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul_si(aux2,aux2,3,prec);
  arb_div_si(aux2,aux2,8,prec);
  arb_add(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
   
  compute_bound_alpha_eps_0(aux2, eps, M, params, prec);
  arb_mul(out,out,aux2,prec);
  arb_div(out,out,twopi,prec);
  arb_mul_si(out,out,4,prec);

  arb_clear(c_ep);
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(lambda);
  arb_clear(twopi);
  return 1;
}


int compute_bound_delta_eps_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  
  arb_init(aux);
  
  compute_c_u0_hat_plus(out, eps, M, params, prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(out,out,aux,prec);
  arb_add_si(out,out,1,prec);
  compute_c_F0(aux, eps, M, params, prec);
  arb_mul(out,out,aux,prec);
  
  arb_clear(aux);
  return 1;
}


int compute_bound_delta_eps_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_c_u0_hat_plus(out, eps, M, params, prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(out,out,aux,prec);
  arb_add_si(out,out,1,prec);
  compute_c_F1(aux, eps, M, params, prec);
  arb_mul(out,out,aux,prec);
  
  compute_cp_p0_plus(aux, eps, M, params, prec);
  arb_pow(aux2,eps,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(aux,aux,lambda,prec);
  compute_bound_delta_eps_0(aux2, eps, M, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_div_si(aux,aux,2,prec);
  arb_div(aux,aux,lambda,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}


int compute_bound_delta_eps_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  
  compute_c_F2(out, eps, M, params, prec);
  
  compute_cpp_p0_zero(aux, eps, M, params, prec);
  arb_pow(aux2,eps,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_add(aux,aux,lambda,prec);
  compute_bound_delta_eps_0(aux2, eps, M, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_div_si(aux,aux,4,prec);
  arb_div(aux,aux,lambda,prec);
  arb_add(out,out,aux,prec);
  
  compute_cp_p0_zero(aux, eps, M, params, prec);
  arb_pow(aux2,eps,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(aux,aux,lambda,prec);
  compute_bound_delta_eps_1(aux2, eps, M, params, prec);
  arb_mul(aux,aux,aux2,prec);
  arb_div_si(aux,aux,2,prec);
  arb_div(aux,aux,lambda,prec);
  arb_add(out,out,aux,prec);
  
  compute_c_u0_hat_plus(aux, eps, M, params, prec);
  arb_pow(aux2,eps,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add_si(aux,aux,1,prec);
  arb_mul(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(lambda);
  return 1;
}


int compute_cT0_1(arb_ptr out, const arb_t inp, void* params, slong prec){
  u_params_struct* ups = params;
  
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
  arb_mul(aux2,inp,inp,prec);
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


int compute_cT0_2(arb_ptr out, const arb_t inp, void* params, slong prec){
  u_params_struct* ups = params;
  
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
  arb_mul(aux2,inp,inp,prec);
  arb_sub(aux2,aux,aux2,prec);
  arb_set_str(argument,"2.5",prec);
  arb_zeta(aux,argument,prec);
  arb_div(aux,aux,aux2,prec);
  arb_sqrt(aux2,pi,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux2,pi,pi,prec);
  arb_mul(aux3,inp,inp,prec);
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


int compute_cT1_1(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  arb_mul_si(aux,inp,2,prec);
  arb_set_str(aux2,"2.5",prec);
  compute_E_x_Cz(out,aux,aux2,params,prec);
  arb_div(out,out,inp,prec);
  
  compute_cp_p0_minus(aux,inp,M,params,prec);
  arb_pow(aux2,inp,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_div(aux,aux,lambda,prec);
  arb_add_si(aux,aux,1,prec);
  arb_mul(out,out,aux,prec);
  
  arb_sqrt_ui(aux,2,prec);
  arb_sub_si(aux,aux,2,prec);
  arb_div_si(aux,aux,8,prec);
  arb_mul(out,out,aux,prec);
  arb_neg(out,out);
  
  arb_mul_si(aux,inp,2,prec);
  arb_set_str(aux2,"1.5",prec);
  compute_E_x_Sz(aux,aux,aux2,params,prec);
  
  arb_sqrt_ui(aux2,2,prec);
  arb_sub_si(aux2,aux2,1,prec);
  arb_div_si(aux2,aux2,2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"1.5",prec);
  compute_E_x_Sz(aux,inp,aux,params,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"-0.5",prec);
  arb_sub(aux,aux,ups->p0,prec);
  arb_pow(aux,inp,aux,prec);
  arb_mul(out,out,aux,prec);

  compute_cp_p0_minus(aux,inp,M,params,prec);
  arb_div(out,out,aux,prec);
  arb_div(out,out,pi,prec);
  
  arb_sqrt_ui(aux,2,prec);
  arb_neg(aux,aux);
  arb_add_si(aux,aux,2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_sqrt(aux2,pi,prec);
  arb_div(aux,aux,aux2,prec);
  arb_div_si(aux,aux,3,prec);
  
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}


int compute_cT1_2(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t pi;
  arb_init(pi);
  arb_const_pi(pi,prec);
  
  arb_set_str(out,"2",prec);
  arb_sqrt(out,out,prec);
  arb_sub_si(out,out,2,prec);
  arb_mul_si(out,out,8,prec);
  arb_div_si(out,out,3,prec);
  arb_div(out,out,pi,prec);
  arb_neg(out,out);
  
  arb_clear(pi);
  return 1;
}


int compute_cT1_3(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  for (m = 1; m <= TANH_MAX; m++){
    arb_set_si(aux,m);
    arb_tanh(aux,aux,prec);
    arb_sqrt(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_set_si(aux2,m);
    arb_sqrt(aux2,aux2,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  arb_set_str(aux,"0.5",prec);
  compute_tail_error(aux,aux,TANH_MAX+1,prec);
  arb_add(out,out,aux,prec);

  arb_mul_si(out,out,5,prec);
  arb_div_si(out,out,4,prec);
  arb_div(out,out,pi,prec);
  
  compute_cp_eps(aux, inp, M, params, prec);
  
  arb_add(out,out,aux,prec);  
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  return 1;
}


int compute_cT1_4(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  arb_t lambda;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
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
  
  compute_cp_p0_plus(aux, inp, M, params, prec);
  arb_pow(aux2,inp,ups->p0,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_div(aux,aux,lambda,prec);
  arb_set_str(aux2,"2.5",prec);
  arb_add(aux,aux,aux2,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,4,prec);
  arb_div(out,out,pi,prec);
  arb_neg(out,out);
  
  arb_sqrt_ui(aux,32,prec);
  arb_sub_si(aux,aux,5,prec);
  arb_set_str(aux2,"-1.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul(aux,aux,inp,prec);
  arb_mul(aux,aux,inp,prec);
  arb_div_si(aux,aux,6,prec);
  arb_div(aux,aux,pi,prec);
  arb_sub(out,out,aux,prec);
  
  arb_sqrt_ui(aux,2,prec);
  arb_sub_si(aux,aux,2,prec);
  arb_set_str(aux2,"0.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,3,prec);
  arb_div_si(aux,aux,4,prec);
  arb_div(aux,aux,pi,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}



int compute_cp_eps(arb_ptr out, const arb_t inp, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  arb_t lambda;
    
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_init(lambda);
  
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_cp_p0_plus(out,inp,M,params,prec);
  arb_div(out,out,lambda,prec);
  arb_pow(aux,inp,ups->p0,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(aux,"2.5",prec);
  arb_add(out,out,aux,prec);
  
  
  arb_mul(aux,pi,pi,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_mul(aux2,inp,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  
  arb_mul(aux,pi,pi,prec);
  arb_mul(aux2,inp,inp,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_div(out,out,aux,prec);
  
  arb_set_str(aux,"1.25",prec);
  arb_zeta(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  
  arb_mul(out,out,inp,prec);
  
  arb_div_si(aux,pi,2,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  
  arb_div_si(out,out,2,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}

int compute_c_p0_SE_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2, aux3;
  arb_t pi;
  arb_t lambda;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
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
  
  arb_mul(out,out,lambda,prec);
  
  arb_sub_si(aux,ups->p0,2,prec);
  arb_neg(aux,aux);
  arb_pow(aux,eps,aux,prec);
  arb_mul(out,out,aux,prec);
  
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
  
  compute_cp_p0_plus(aux2, eps, M, params, prec);
  arb_mul_si(aux2,aux2,2,prec);
  compute_cpp_p0_minus(aux3, eps, M, params, prec);
  arb_add(aux2,aux2,aux3,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_div(aux,pi,eps,prec);
  arb_log(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  
  arb_sqrt(aux,eps,prec);
  arb_mul(out,out,aux,prec);
  
  arb_div_si(out,out,2,prec);
  arb_div(out,out,pi,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}


int compute_c_p0_SE_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  arb_t lambda;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
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
  
  
  compute_cp_p0_plus(aux, eps, M, params, prec);
  arb_mul_si(aux,aux,2,prec);
  compute_cpp_p0_plus(aux2, eps, M, params, prec);
  arb_add(aux,aux,aux2,prec);
  arb_div(aux,aux,lambda,prec);
  arb_set_str(aux2,"19",prec);
  arb_div_si(aux2,aux2,8,prec);
  arb_add(aux,aux,aux2,prec);
  arb_mul(out,out,aux,prec);
  
  arb_div_si(out,out,4,prec);
  arb_div(out,out,pi,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}


int compute_c_p0_T2_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){  
  arb_t aux, aux2, aux3;
  arb_t pi;
  arb_t lambda;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  // Third line
  arb_mul_si(out,pi,2,prec);
  arb_set_si(aux,3);
  arb_sqrt(aux,aux,prec);
  arb_div(out,out,aux,prec);
  arb_div(out,out,eps,prec);
  arb_log(out,out,prec);
  arb_div_si(out,out,2,prec);
  arb_div(out,out,pi,prec);
  arb_div(out,out,pi,prec);
  arb_set_si(aux,83);
  arb_div_si(aux,aux,216,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_add(out,out,aux,prec);
  arb_mul(out,out,eps,prec);
  arb_mul(out,out,eps,prec);
  
  arb_set_si(aux,4);
  arb_div_si(aux,aux,3,prec);
  arb_log(aux,aux,prec);
  arb_add(out,out,aux,prec);
  arb_set_si(aux,17);
  arb_div_si(aux,aux,9,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"2.5",prec);
  arb_zeta(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div(out,out,lambda,prec);
  arb_div_si(out,out,4,prec);
  arb_div(out,out,pi,prec);
  
  arb_div(aux,pi,eps,prec);
  arb_log(aux,aux,prec);
  arb_set_str(aux2,"0.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_abs(aux2,aux2);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  compute_cpp_p0_minus(aux,eps,M,params,prec);
  arb_mul(out,out,aux,prec);
  
  // middle line
  arb_div(aux,pi,eps,prec);
  arb_log(aux,aux,prec);
  arb_mul(aux,aux,eps,prec);
  arb_mul(aux,aux,eps,prec);
  arb_mul(aux2,pi,pi,prec);
  arb_mul_si(aux2,aux2,3,prec);
  arb_div_si(aux2,aux2,2,prec);
  arb_add(aux,aux,aux2,prec);
  arb_set_str(aux2,"-1.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_abs(aux2,aux2);
  arb_mul(aux,aux,aux2,prec);
  arb_div_si(aux,aux,3,prec);
  
  arb_div(aux2,pi,eps,prec);
  arb_log(aux2,aux2,prec);
  arb_set_str(aux3,"0.5",prec);
  arb_zeta(aux3,aux3,prec);
  arb_abs(aux3,aux3);
  arb_mul(aux2,aux2,aux3,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_add(aux,aux,aux2,prec);
  
  // first line
  arb_set_si(aux2,4);
  arb_div_si(aux2,aux2,3,prec);
  arb_log(aux2,aux2,prec);
  arb_set_si(aux3,16);
  arb_div_si(aux3,aux3,3,prec);
  arb_add(aux2,aux2,aux3,prec);
  arb_mul(aux2,aux2,eps,prec);
  arb_mul(aux2,aux2,eps,prec);
  
  arb_set_si(aux3,4);
  arb_div_si(aux3,aux3,3,prec);
  arb_log(aux3,aux3,prec);
  arb_mul_si(aux3,aux3,12,prec);
  arb_sub_si(aux3,aux3,1,prec);
  arb_mul(aux3,aux3,pi,prec);
  arb_mul(aux3,aux3,pi,prec);
  arb_div_si(aux3,aux3,3,prec);
  arb_add(aux2,aux2,aux3,prec);
  
  arb_set_str(aux3,"4.5",prec);
  arb_zeta(aux3,aux3,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_div(aux2,aux2,lambda,prec);
  arb_div_si(aux2,aux2,16,prec);
  arb_div(aux2,aux2,pi,prec);
  arb_div(aux2,aux2,pi,prec);
  arb_div(aux2,aux2,pi,prec);
  arb_add(aux,aux,aux2,prec);
  
  
  compute_cp_p0_minus(aux2,eps,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_sqrt(aux,eps,prec);
  arb_mul(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}


int compute_c_p0_T2_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  
  arb_const_pi(pi,prec);
  
  arb_sqrt_ui(out,50,prec);
  arb_neg(out,out);
  arb_add_si(out,out,2,prec);
  arb_mul_si(out,out,2,prec);
  arb_div_si(out,out,3,prec);
  arb_add(out,out,pi,prec);
  arb_set_str(aux,"1",prec);
  arb_asinh(aux,aux,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  arb_div_si(out,out,3,prec);
  compute_cpp_p0_plus(aux,eps,M,params,prec);
  arb_mul(out,out,aux,prec);
  
  arb_set_str(aux,"1",prec);
  arb_asinh(aux,aux,prec);
  arb_sqrt_ui(aux2,2,prec);
  arb_sub(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(aux,aux,pi,prec);
  compute_cp_p0_minus(aux2,eps,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  arb_mul_si(out,out,4,prec);
  arb_div(out,out,pi,prec);
  
  compute_c_p0_SE_1(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  
  compute_c_p0_T2_0(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  return 1;
}


int compute_c_p0_T2_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  
  arb_const_pi(pi,prec);
  
  arb_set_str(out,"2",prec);
  arb_log(out,out,prec);
  arb_mul_si(out,out,-2,prec);
  arb_add_si(out,out,3,prec);
  
  compute_cpp_p0_plus(aux,eps,M,params,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,3,prec);
  
  arb_set_str(aux,"2",prec);
  arb_log(aux,aux,prec);
  arb_neg(aux,aux);
  arb_add_si(aux,aux,1,prec);
  
  compute_cp_p0_plus(aux2,eps,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,2,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul_si(out,out,8,prec);
  arb_div(out,out,pi,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  return 1;
}


int compute_cT2_1(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  arb_t lambda;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  arb_mul_si(out,pi,2,prec);
  arb_div(out,out,eps,prec);
  arb_set_si(aux,3);
  arb_sqrt(aux,aux,prec);
  arb_div(out,out,aux,prec);
  arb_log(out,out,prec);
  
  arb_set_si(aux,83);
  arb_div_si(aux,aux,108,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"2.5",prec);
  arb_zeta(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_mul_si(out,out,3,prec);
  arb_div_si(out,out,2,prec);
  
  arb_set_si(aux,4);
  arb_div_si(aux,aux,3,prec);
  arb_log(aux,aux,prec);
  arb_set_si(aux2,16);
  arb_div_si(aux2,aux2,3,prec);
  arb_add(aux,aux,aux2,prec);
  arb_set_si(aux2,9);
  arb_div_si(aux2,aux2,2,prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul(out,out,eps,prec);
  arb_div(out,out,lambda,prec);
  arb_div_si(out,out,16,prec);
  arb_div(out,out,pi,prec);
  arb_div(out,out,pi,prec);
  arb_div(out,out,pi,prec);
  
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}

int compute_cT2_2(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux;
  arb_init(aux);
  
  compute_c_p0_SE_2(out,eps,M,params,prec);
  compute_cT2_2_0(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  arb_mul(out,out,eps,prec);
  
  arb_clear(aux);
  return 1;
}


int compute_cT2_2_0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2, aux3;
  arb_t pi;
  arb_t lambda;
  
  int m;
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(aux3);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  compute_cpp_p0_plus(out, eps, M, params, prec);
  arb_pow(aux,eps,ups->p0,prec);
  arb_mul(out,out,aux,prec);
  arb_div(out,out,lambda,prec);
  arb_set_str(aux,"0.75",prec);
  arb_add(out,out,aux,prec);
  arb_set_str(aux,"2.5",prec);
  arb_zeta(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_set_str(aux,"2",prec);
  arb_div(aux,aux,pi,prec);
  arb_sqrt(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_mul_si(out,out,2,prec);
  arb_div(out,out,pi,prec);
  arb_div(out,out,pi,prec);
  
  arb_mul(aux,eps,eps,prec);
  arb_mul_si(aux,aux,9,prec);
  arb_div_si(aux,aux,2,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_set_str(aux2,"-1.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_div_si(aux2,aux2,4,prec);
  arb_sub(aux,aux,aux2,prec);
  compute_cp_p0_plus(aux2,eps,M,params,prec);
  arb_pow(aux3,eps,ups->p0,prec);
  arb_mul(aux2,aux2,aux3,prec);
  arb_div(aux2,aux2,lambda,prec);
  arb_mul_si(aux2,aux2,2,prec);
  arb_add_si(aux2,aux2,1,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_str(aux2,"4.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_set_str(aux2,"2",prec);
  arb_div(aux2,aux2,pi,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_div_si(aux,aux,4,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"-1.5",prec);
  arb_zeta(aux,aux,prec);
  arb_div_si(aux,aux,4,prec);
  arb_sub(out,out,aux,prec);
  
  arb_set_str(aux,"0.5",prec);
  arb_zeta(aux,aux,prec);
  arb_abs(aux,aux);
  arb_set_str(aux2,"2",prec);
  arb_div(aux2,aux2,pi,prec);
  arb_sqrt(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_mul_si(aux,aux,4,prec);
  arb_div_si(aux,aux,9,prec);
  arb_div(aux,aux,pi,prec);
  arb_div(aux,aux,pi,prec);
  arb_add(out,out,aux,prec);
  
  arb_div_si(out,out,2,prec);
  arb_div(out,out,pi,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(aux3);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}


int compute_cT2_3(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  arb_t pi;
  arb_t lambda;
  
  
  arb_init(aux);
  arb_init(aux2);
  arb_init(pi);
  arb_init(lambda);
  arb_set_str(out,"0",prec);
  arb_const_pi(pi,prec);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
    
  arb_set_si(out,4);
  arb_div_si(out,out,3,prec);
  arb_log(out,out,prec);
  arb_mul_si(out,out,12,prec);
  arb_sub_si(out,out,1,prec);
  
  arb_set_str(aux,"4.5",prec);
  arb_zeta(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_div_si(out,out,9,prec);
  
  arb_set_si(aux,4);
  arb_div_si(aux,aux,3,prec);
  arb_log(aux,aux,prec);
  arb_set_si(aux2,17);
  arb_div_si(aux2,aux2,9,prec);
  arb_add(aux,aux,aux2,prec);
  arb_set_str(aux2,"2.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  arb_div(out,out,lambda,prec);
  arb_div_si(out,out,-2,prec);
  
  arb_div(aux,pi,eps,prec);
  arb_log(aux,aux,prec);
  arb_set_str(aux2,"0.5",prec);
  arb_zeta(aux2,aux2,prec);
  arb_abs(aux2,aux2);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  
  arb_mul_si(out,out,3,prec);
  arb_div_si(out,out,8,prec);
  arb_div(out,out,pi,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  arb_clear(pi);
  arb_clear(lambda);
  return 1;
}






// (XXX) in the paper. We prove c'_{p0,0}*eps^p0 - lambda/2 < 0
int compute_validation_du0(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t c_ep;
  arb_t aux;
  arb_t lambda;
  
  
  arb_init(c_ep);
  arb_init(aux);
  arb_init(lambda);
  
  // initialize lambda
  arb_const_pi(lambda,prec);
  arb_div_si(lambda,lambda,8,prec);
  arb_sqrt(lambda,lambda,prec);
  
  arb_pow(aux,eps,ups->p0,prec);
  compute_cp_p0_zero(c_ep,eps,M,params,prec);
  arb_mul(out,aux,c_ep,prec);
  arb_div_si(aux,lambda,2,prec);
  arb_sub(out,out,aux,prec);
  
  arb_clear(c_ep);
  arb_clear(aux);
  arb_clear(lambda);
  return 1;
}


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




// (XXX) in the paper. 
// We prove (c_{p0} + c'_{p0,-}c_{T_1}^{1} + c'_{p0,+}c_{T_1}^{2}))\ep^{p0-1/2} + c_{T_1}^{3} \ep - c_{T_1}^{4} < 0
int compute_validation_T1(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  
  arb_init(aux);
  arb_init(aux2);
  
  compute_cT1_1(out,eps,M,params,prec);
  compute_cp_p0_minus(aux,eps,M,params,prec);
  arb_mul(out,out,aux,prec);
  compute_cT1_2(aux,eps,M,params,prec);
  compute_cp_p0_plus(aux2,eps,M,params,prec);
  arb_mul(aux,aux,aux2,prec);
  arb_add(out,out,aux,prec);
  compute_c_p0(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"0.5",prec);
  arb_sub(aux,ups->p0,aux,prec);
  arb_pow(aux,eps,aux,prec);
  arb_mul(out,out,aux,prec);
  
  compute_cT1_3(aux,eps,M,params,prec);
  arb_mul(aux,aux,eps,prec);
  arb_add(out,out,aux,prec);
  
  compute_cT1_4(aux,eps,M,params,prec);
  arb_sub(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  return 1;
}


// (A.11) in the paper. 
// We prove XXX (FILL THIS)
// TODO(Javi): untested
int compute_validation_T2(arb_ptr out, const arb_t eps, int M, void* params, slong prec){
  u_params_struct* ups = params;
  
  arb_t aux, aux2;
  
  arb_init(aux);
  arb_init(aux2);
  
  compute_c_p0(out,eps,M,params,prec);
  compute_c_p0_T2_1(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  compute_c_p0_T2_2(aux,eps,M,params,prec);
  arb_add(out,out,aux,prec);
  
  arb_set_str(aux,"0.5",prec);
  arb_sub(aux,ups->p0,aux,prec);
  arb_pow(aux,eps,aux,prec);
  arb_mul(out,out,aux,prec);
  
  compute_cT2_1(aux,eps,M,params,prec);
  compute_cT2_2(aux2,eps,M,params,prec);
  arb_add(aux,aux,aux2,prec);
  arb_mul(aux,aux,eps,prec);
  arb_add(out,out,aux,prec);
    
  compute_cT2_3(aux,eps,M,params,prec);
  arb_sub(out,out,aux,prec);
  
  arb_clear(aux);
  arb_clear(aux2);
  return 1;
}


int compute_nu(arb_ptr out, void* params, slong prec){
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
  
  // We first calculate 4*Lu_0(0)
  
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
  arb_set_str(aux,"0.0000359",prec);   // Hardcode eps_0
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

int compute_mu(arb_ptr out, void* params, slong prec){
  compute_nu(out,params,prec);
  arb_neg(out,out);
  arb_add_si(out,out,1,prec);
  arb_sqrt(out,out,prec);
  arb_neg(out,out);
  arb_add_si(out,out,1,prec);
  return 1;
}
