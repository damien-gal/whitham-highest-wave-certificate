/* basic_functions.c: shared evaluators of the computer-assisted proof.
 * See basic_functions.h for the list of routines. The Fourier sums are
 * truncated after n = TANH_MAX and their tails are bounded by the Fourier-tail
 * lemma of Appendix B (compute_tail_error).
 */
#include "basic_functions.h"
#include "assert.h"

void arb_clausen_sin_direct(arb_t out, const arb_t s, const arb_t x, slong prec){
  if(clausen_cached(out,s,x,1,prec)) return;
  if(!arb_is_int(s)) { arb_indeterminate(out); return; }

  arb_t cx,sx;
  acb_t poly,polar;
  acb_t complex_s;
  
  arb_init(cx);
  arb_init(sx);
  acb_init(poly);
  acb_init(polar);
  acb_init(complex_s);

  acb_set_arb(complex_s,s);  
  arb_sin_cos(sx,cx,x,prec);
  acb_set_arb_arb(polar,cx,sx);
  acb_polylog(poly,complex_s,polar,prec);
  acb_get_imag(out,poly);
  
  arb_clear(cx);
  arb_clear(sx);
  acb_clear(poly);
  acb_clear(polar);
  acb_clear(complex_s);
  
}

void arb_clausen_cos(arb_t out, const arb_t s, const arb_t x, slong prec){
  if(clausen_cached(out,s,x,0,prec)) return;
  if(!arb_is_int(s)) { arb_indeterminate(out); return; }

  arb_t clx,slx;
  arb_t crx,srx;
  arb_t outl,outr,outz;
  acb_t polyl,polarl;
  acb_t polyr,polarr;
  acb_t polyz,polarz;
  acb_t complex_s;
  
  arf_t left,right,zero;
  arb_t minus_one;
  arb_t pi;
  
  int close_to_pi = 0; // 1 if x is close to pi, 0 otherwise
  
  arf_init(left); arf_init(right); arf_init(zero); arb_init(minus_one); arb_init(pi);
  
  arf_zero(zero); 
  arb_one(minus_one); arb_neg(minus_one,minus_one);
  arb_const_pi(pi,prec);
  
  arb_init(clx);
  arb_init(slx);
  acb_init(polyl);
  acb_init(polarl);
  acb_init(complex_s);
  
  arb_init(crx);
  arb_init(srx);
  acb_init(polyr);
  acb_init(polarr);
  
  acb_init(polyz);
  acb_init(polarz);
  
  arb_init(outl); arb_init(outr); arb_init(outz);
  
  arb_get_lbound_arf(left,x,prec);
  arb_get_ubound_arf(right,x,prec);
  
  // Checks with endpoints
  
  // We take out negative values
  arf_max(left,zero,left);
  
  // We check if we are close to pi. If so, the final result is the union of what we got with the value at pi.
  if (arb_overlaps(x,pi)){
    close_to_pi = 1;
  }

  acb_set_arb(complex_s,s); 
  
  arb_sin_cos_arf_generic(slx,clx,left,prec);
  acb_set_arb_arb(polarl,clx,slx);
  acb_polylog(polyl,complex_s,polarl,prec);
  acb_get_real(outl,polyl);
  
  arb_sin_cos_arf_generic(srx,crx,right,prec);
  acb_set_arb_arb(polarr,crx,srx);
  acb_polylog(polyr,complex_s,polarr,prec);
  acb_get_real(outr,polyr);
  
  arb_union(out,outl,outr,prec);
  
  if (close_to_pi){
    acb_set_arb(polarz,minus_one);
    acb_polylog(polyz,complex_s,polarz,prec);
    acb_get_real(outz,polyz);
    arb_union(out,out,outz,prec);
  }
  
  arf_clear(left); arf_clear(right); arf_clear(zero); arb_clear(minus_one); arb_clear(pi);
  
  arb_clear(outl); arb_clear(outr); arb_clear(outz);
  
  arb_clear(clx);
  arb_clear(slx);
  acb_clear(polyl);
  acb_clear(polarl);
  acb_clear(complex_s);
  
  acb_clear(polyz);
  acb_clear(polarz);
  
  arb_clear(crx);
  arb_clear(srx);
  acb_clear(polyr);
  acb_clear(polarr);
  
}


void arb_clausen_sin(arb_t out, const arb_t s, const arb_t x, slong prec){
  // Implemented through arb_clausen_cos using its monotonicity
  
  arb_t mid;
  arb_t fmid;
  arb_t deriv;
  arb_t sderiv;
  
  arb_init(mid);
  arb_init(fmid);
  arb_init(deriv);
  arb_init(sderiv);
  
  arb_get_mid_arb(mid,x);
  
  arb_clausen_sin_direct(fmid,s,mid,prec);
  arb_sub_si(sderiv,s,1,prec);
  arb_clausen_cos(deriv,sderiv,x,prec);
  
  arb_sub(out,x,mid,prec);
  arb_mul(out,out,deriv,prec);
  arb_add(out,out,fmid,prec);
  
  arb_clear(mid);
  arb_clear(fmid);
  arb_clear(deriv);
  arb_clear(sderiv);
  
}

// The Fourier-tail lemma of Appendix B with p = exponent and n_0 = N0:
//   sum_{n>=n_0} n^p (1 - sqrt(tanh n)) <= 2 n_0^{min(p,0)} / (e^{n_0 q} - e^{(n_0-1) q}),
//   q = 2 - max(p,0) log(n_0)/n_0 > 0.
// Evaluated with ball arithmetic in the exponent, so the result bounds the
// tail for every exponent in the ball. The positivity of q is asserted.
void compute_tail_error(arb_ptr out, const arb_t exponent, const int N0, slong prec){

  arb_t minp,maxp;
  arb_t pN0;
  arb_t aux;
  arb_t zero;
  
  
  arb_init(minp);
  arb_init(maxp);
  arb_init(pN0);
  arb_init(aux);
  arb_init(zero);
  
  arb_zero(zero);
  
  // min(p,0) and max(p,0) as balls
  arb_min(minp,zero,exponent,prec);
  arb_max(maxp,zero,exponent,prec);
  
  arb_log_ui(pN0,N0,prec);
  arb_div_si(pN0,pN0,N0,prec);
  arb_mul(pN0,pN0,maxp,prec);
  arb_neg(pN0,pN0);
  arb_add_si(pN0,pN0,2,prec);
  
  assert(arb_is_positive(pN0));
  
  arb_set_si(aux,N0);
  arb_pow(out,aux,minp,prec);
  arb_mul_si(aux,pN0,-N0,prec);
  arb_exp(aux,aux,prec);
  arb_mul(out,out,aux,prec);
  arb_neg(aux,pN0);
  arb_exp(aux,aux,prec);
  arb_neg(aux,aux);
  arb_add_si(aux,aux,1,prec);
  arb_div(out,out,aux,prec);
  arb_mul_si(out,out,2,prec);
  
  arb_clear(minp);
  arb_clear(maxp);
  arb_clear(pN0);
  arb_clear(aux);
  arb_clear(zero);
}


// Populates a u_params_struct with the right parameters for u0
void populate_parameters(u_params_struct* up, slong prec){
  up->N0 = 6; up->N1 = 2; up->N2 = 11;

  arb_init(up->p0); arb_init(up->p1);

  arb_set_str(up->p0,"0.61120159320355359",prec);
  arb_set_str(up->p1,"2.7624012994582795",prec);
  
  up->a0k = _arb_vec_init(up->N0);
  up->a1k = _arb_vec_init(up->N1);
  up->bi = _arb_vec_init(up->N2);

  arb_set_str(up->a0k,"0.25",prec);
  arb_set_str(up->a0k+1,"-0.09186529326005956",prec);
  arb_set_str(up->a0k+2,"-0.004027567598029841",prec);
  arb_set_str(up->a0k+3,"0.002449055302062304",prec);
  arb_set_str(up->a0k+4,"0.005300258772850143",prec);
  arb_set_str(up->a0k+5,"0.0005045620668399171",prec);
  
  arb_set_str(up->a1k,"-0.015079514055110828",prec);
  arb_set_str(up->a1k+1,"-0.0015036862476077348",prec);
  
  arb_set_str(up->bi,"0.008560417336473005",prec);
  arb_set_str(up->bi+1,"0.0006326823377277547",prec);
  arb_set_str(up->bi+2,"0.00004973972422133443",prec);
  arb_set_str(up->bi+3,"4.886926723089561e-6",prec);
  arb_set_str(up->bi+4,"5.532358724209955e-7",prec);
  arb_set_str(up->bi+5,"5.719481296653965e-8",prec);
  arb_set_str(up->bi+6,"1.4354134439015963e-8",prec);
  arb_set_str(up->bi+7,"-9.97622130913055e-9",prec);
  arb_set_str(up->bi+8,"5.828491310021984e-9",prec);
  arb_set_str(up->bi+9,"-1.847883895374688e-8",prec);
  arb_set_str(up->bi+10,"-4.463362465423652e-8",prec);
}

int compute_u0(arb_ptr out, const arb_t inp, void* params, slong prec){
  u_params_struct* up = params;
  arb_t a0kterm,a1kterm,biterm;
  arb_t aux;
  arb_t clausen;
  arb_t zeta;
  arb_t argument;
  
  int i;
  
  arb_init(a0kterm); arb_init(a1kterm); arb_init(biterm);
  arb_zero(a0kterm); arb_zero(a1kterm); arb_zero(biterm);

  arb_init(aux); arb_init(clausen); arb_init(zeta); arb_init(argument);
  
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_zeta(zeta,argument,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_sub(aux,zeta,clausen,prec);
    arb_mul(aux,aux,up->a0k+i,prec);
    arb_add(a0kterm,a0kterm,aux,prec);
  }

  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"1.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_zeta(zeta,argument,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_sub(aux,zeta,clausen,prec);
    arb_mul(aux,aux,up->a1k+i,prec);
    arb_add(a1kterm,a1kterm,aux,prec);
  }
  
  for (i = 1; i <= up->N2; i++){
    arb_mul_si(aux,inp,i,prec);
    arb_cos(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    arb_mul(aux,aux,up->bi+i-1,prec);
    arb_add(biterm,biterm,aux,prec);
  }
  
  arb_add(out,a0kterm,a1kterm,prec);
  arb_add(out,out,biterm,prec);
  
  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);
  
  arb_clear(aux);
  arb_clear(clausen);
  arb_clear(zeta);
  arb_clear(argument);
  return 0;
}

int compute_du0(arb_ptr out, const arb_t inp, void* params, slong prec){  
  u_params_struct* up = params;
  arb_t a0kterm,a1kterm,biterm;
  arb_t aux;
  arb_t clausen;
  arb_t argument;
  
  int i;
  
  arb_init(a0kterm); arb_init(a1kterm); arb_init(biterm);
  arb_zero(a0kterm); arb_zero(a1kterm); arb_zero(biterm);

  arb_init(aux); arb_init(clausen); arb_init(argument);
  
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_clausen_sin(clausen,argument,inp,prec);
    arb_mul(clausen,clausen,up->a0k+i,prec);
    arb_add(a0kterm,a0kterm,clausen,prec);
  }
  
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_clausen_sin(clausen,argument,inp,prec);
    arb_mul(clausen,clausen,up->a1k+i,prec);
    arb_add(a1kterm,a1kterm,clausen,prec);
  }
  
  for (i = 1; i <= up->N2; i++){
    arb_mul_si(aux,inp,i,prec);
    arb_sin(aux,aux,prec);
    arb_mul_si(aux,aux,i,prec);
    arb_mul(aux,aux,up->bi+i-1,prec);
    arb_neg(aux,aux);
    arb_add(biterm,biterm,aux,prec);
  }
  
  
  arb_add(out,a0kterm,a1kterm,prec);
  arb_add(out,out,biterm,prec);
  
  
  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);
  
  arb_clear(aux);
  arb_clear(clausen);
  arb_clear(argument);
  return 0;
}

int compute_ddu0(arb_ptr out, const arb_t inp, void* params, slong prec){  
  u_params_struct* up = params;
  arb_t a0kterm,a1kterm,biterm;
  arb_t aux;
  arb_t clausen;
  arb_t argument;
  
  int i;
  
  arb_init(a0kterm); arb_init(a1kterm); arb_init(biterm);
  arb_zero(a0kterm); arb_zero(a1kterm); arb_zero(biterm);
  
  arb_init(aux); arb_init(clausen); arb_init(argument);
  
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_mul(clausen,clausen,up->a0k+i,prec);
    arb_add(a0kterm,a0kterm,clausen,prec);
  }
  
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"-0.5",prec);
    arb_add(argument,argument,aux,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_mul(clausen,clausen,up->a1k+i,prec);
    arb_add(a1kterm,a1kterm,clausen,prec);
  }
  
  for (i = 1; i <= up->N2; i++){
    arb_mul_si(aux,inp,i,prec);
    arb_cos(aux,aux,prec);
    arb_mul_si(aux,aux,i,prec);
    arb_mul_si(aux,aux,i,prec);
    arb_mul(aux,aux,up->bi+i-1,prec);
    arb_neg(aux,aux);
    arb_add(biterm,biterm,aux,prec);
  }
  
  
  arb_add(out,a0kterm,a1kterm,prec);
  arb_add(out,out,biterm,prec);
  
  
  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);
  
  arb_clear(aux);
  arb_clear(clausen);
  arb_clear(argument);
  return 0;
}


int compute_Lu0(arb_ptr out, const arb_t inp, void* params, slong prec){
  u_params_struct* up = params;
  arb_t a0kterm,a1kterm,biterm;
  arb_t aux, aux2;
  arb_t clausen;
  arb_t zeta;
  arb_t argument;
  arb_t error_tail;
  arb_t pm_one;
  
  int i,m;
  
  arb_init(a0kterm); arb_init(a1kterm); arb_init(biterm);
  arb_zero(a0kterm); arb_zero(a1kterm); arb_zero(biterm);
  
  arb_init(aux); arb_init(aux2);
  arb_init(clausen); arb_init(zeta); arb_init(argument); arb_init(error_tail); arb_init(pm_one);
  arb_zero_pm_one(pm_one);
  
  // Main term
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"2",prec);
    arb_add(argument,argument,aux,prec);
    arb_zeta(zeta,argument,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_sub(aux,zeta,clausen,prec);
    arb_mul(aux,aux,up->a0k+i,prec);
    arb_add(a0kterm,a0kterm,aux,prec);
  }
  
  // Main term
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"2",prec);
    arb_add(argument,argument,aux,prec);
    arb_zeta(zeta,argument,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_sub(aux,zeta,clausen,prec);
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
    arb_mul_si(error_tail,error_tail,2,prec);
    arb_mul(error_tail,error_tail,up->a0k+i,prec);
    
    for (m = 1; m <= TANH_MAX; m++){
      
      arb_set_si(aux2,m);
      arb_tanh(aux2,aux2,prec);
      arb_sqrt(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      
      arb_mul_si(aux,inp,m,prec);
      arb_cos(aux,aux,prec);
      arb_sub_si(aux,aux,1,prec);
      arb_neg(aux,aux);
      
      arb_mul(aux,aux,aux2,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a0k+i,prec);
      arb_add(a0kterm,a0kterm,aux,prec);
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
    arb_mul_si(error_tail,error_tail,2,prec);
    arb_mul(error_tail,error_tail,up->a1k+i,prec);
    
    for (m = 1; m <= TANH_MAX; m++){
      
      arb_set_si(aux2,m);
      arb_tanh(aux2,aux2,prec);
      arb_sqrt(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      
      arb_mul_si(aux,inp,m,prec);
      arb_cos(aux,aux,prec);
      arb_sub_si(aux,aux,1,prec);
      arb_neg(aux,aux);
      
      arb_mul(aux,aux,aux2,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a1k+i,prec);
      arb_add(a1kterm,a1kterm,aux,prec);
    }
    // Tail N1
    arb_mul(error_tail,error_tail,pm_one,prec);
    arb_add(a1kterm,a1kterm,error_tail,prec);
  }
  
  for (i = 1; i <= up->N2; i++){
    // b_k*(cos(k*x)-1)
    arb_mul_si(aux,inp,i,prec);
    arb_cos(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    arb_mul(aux,aux,up->bi+i-1,prec);
    
    //sqrt(tanh(k)/k)
    arb_set_si(aux2,i);
    arb_tanh(aux2,aux2,prec);
    arb_div_si(aux2,aux2,i,prec);
    arb_sqrt(aux2,aux2,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_add(biterm,biterm,aux,prec);
  }
  
  
  arb_add(out,a0kterm,a1kterm,prec);
  arb_add(out,out,biterm,prec);
  
  
  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);
  
  arb_clear(aux); arb_clear(aux2);
  arb_clear(clausen);
  arb_clear(zeta);
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pm_one);
  return 0;
}


int compute_dLu0(arb_ptr out, const arb_t inp, void* params, slong prec){
  u_params_struct* up = params;
  arb_t a0kterm,a1kterm,biterm;
  arb_t aux, aux2;
  arb_t clausen;
  arb_t zeta;
  arb_t argument;
  arb_t error_tail;
  arb_t pm_one;
  
  int i,m;
  
  
  arb_init(a0kterm); arb_init(a1kterm); arb_init(biterm);
  arb_zero(a0kterm); arb_zero(a1kterm); arb_zero(biterm);
  
  arb_init(aux); arb_init(aux2);
  arb_init(clausen); arb_init(zeta); arb_init(argument); arb_init(error_tail); arb_init(pm_one);
  arb_zero_pm_one(pm_one);
  
  // Main term
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"1",prec);
    arb_add(argument,argument,aux,prec);
    arb_clausen_sin(clausen,argument,inp,prec);
    arb_mul(aux,clausen,up->a0k+i,prec);
    arb_add(a0kterm,a0kterm,aux,prec);
  }
  
  // Main term
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"1",prec);
    arb_add(argument,argument,aux,prec);
    arb_clausen_sin(clausen,argument,inp,prec);
    arb_mul(aux,clausen,up->a1k+i,prec);
    arb_add(a1kterm,a1kterm,aux,prec);
  }
  
  // Correction N0
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_set_str(aux,"1",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    
    compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
    arb_mul(error_tail,error_tail,up->a0k+i,prec);
    
    for (m = 1; m <= TANH_MAX; m++){
      
      arb_set_si(aux2,m);
      arb_tanh(aux2,aux2,prec);
      arb_sqrt(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      
      arb_mul_si(aux,inp,m,prec);
      arb_sin(aux,aux,prec);
      
      arb_mul(aux,aux,aux2,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a0k+i,prec);
      arb_add(a0kterm,a0kterm,aux,prec);
    }
    // Tail N0
    arb_mul(error_tail,error_tail,pm_one,prec);
    arb_add(a0kterm,a0kterm,error_tail,prec);
  }
  
  // Correction N1
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_set_str(aux,"1",prec);
    arb_add(argument,argument,aux,prec);
    arb_neg(argument,argument);
    
    compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
    arb_mul(error_tail,error_tail,up->a1k+i,prec);
    for (m = 1; m <= TANH_MAX; m++){
      
      arb_set_si(aux2,m);
      arb_tanh(aux2,aux2,prec);
      arb_sqrt(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      
      arb_mul_si(aux,inp,m,prec);
      arb_sin(aux,aux,prec);
      
      arb_mul(aux,aux,aux2,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a1k+i,prec);
      arb_add(a1kterm,a1kterm,aux,prec);
    }
    // Tail N1
    arb_mul(error_tail,error_tail,pm_one,prec);
    arb_add(a1kterm,a1kterm,error_tail,prec);
  }
  
  for (i = 1; i <= up->N2; i++){
    // b_k*(cos(k*x)-1)
    arb_mul_si(aux,inp,i,prec);
    arb_sin(aux,aux,prec);
    arb_neg(aux,aux);
    arb_mul(aux,aux,up->bi+i-1,prec);
    
    //sqrt(k*tanh(k))
    arb_set_si(aux2,i);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,i,prec);
    arb_sqrt(aux2,aux2,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_add(biterm,biterm,aux,prec);
  }
  
  
  arb_add(out,a0kterm,a1kterm,prec);
  arb_add(out,out,biterm,prec);
  
  
  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);
  
  arb_clear(aux); arb_clear(aux2);
  arb_clear(clausen);
  arb_clear(zeta);
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pm_one);
  return 0;
}


int compute_ddLu0(arb_ptr out, const arb_t inp, void* params, slong prec){
  u_params_struct* up = params;
  arb_t a0kterm,a1kterm,biterm;
  arb_t aux, aux2;
  arb_t clausen;
  arb_t zeta;
  arb_t argument;
  arb_t error_tail;
  arb_t pm_one;
  
  int i,m;
  
  
  arb_init(a0kterm); arb_init(a1kterm); arb_init(biterm);
  arb_zero(a0kterm); arb_zero(a1kterm); arb_zero(biterm);
  
  arb_init(aux); arb_init(aux2);
  arb_init(clausen); arb_init(zeta); arb_init(argument); arb_init(error_tail); arb_init(pm_one);
  arb_zero_pm_one(pm_one);
  
  // Main term
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_mul(aux,clausen,up->a0k+i,prec);
    arb_add(a0kterm,a0kterm,aux,prec);
  }
  
  // Main term
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_clausen_cos(clausen,argument,inp,prec);
    arb_mul(aux,clausen,up->a1k+i,prec);
    arb_add(a1kterm,a1kterm,aux,prec);
  }
  
  // Correction N0
  for (i = 0; i < up->N0; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_neg(argument,argument);
    
    compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
    arb_mul(error_tail,error_tail,up->a0k+i,prec);
    for (m = 1; m <= TANH_MAX; m++){
      
      arb_set_si(aux2,m);
      arb_tanh(aux2,aux2,prec);
      arb_sqrt(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      
      arb_mul_si(aux,inp,m,prec);
      arb_cos(aux,aux,prec);
      
      arb_mul(aux,aux,aux2,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a0k+i,prec);
      arb_add(a0kterm,a0kterm,aux,prec);
    }
    // Tail N0
    arb_mul(error_tail,error_tail,pm_one,prec);
    arb_add(a0kterm,a0kterm,error_tail,prec);
  }
  
  // Correction N1
  for (i = 0; i < up->N1; i++){
    arb_mul_si(argument,up->p0,i,prec);
    arb_add(argument,argument,up->p1,prec);
    arb_neg(argument,argument);
    
    compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
    arb_mul(error_tail,error_tail,up->a1k+i,prec);
    for (m = 1; m <= TANH_MAX; m++){
      
      arb_set_si(aux2,m);
      arb_tanh(aux2,aux2,prec);
      arb_sqrt(aux2,aux2,prec);
      arb_sub_si(aux2,aux2,1,prec);
      
      arb_mul_si(aux,inp,m,prec);
      arb_cos(aux,aux,prec);
      
      arb_mul(aux,aux,aux2,prec);
      
      arb_set_si(aux2,m);
      arb_pow(aux2,aux2,argument,prec);
      
      arb_mul(aux,aux,aux2,prec);
      arb_mul(aux,aux,up->a1k+i,prec);
      arb_add(a1kterm,a1kterm,aux,prec);
    }
    // Tail N1
    arb_mul(error_tail,error_tail,pm_one,prec);
    arb_add(a1kterm,a1kterm,error_tail,prec);
  }
  
  for (i = 1; i <= up->N2; i++){
    // b_k*(cos(k*x)-1)
    arb_mul_si(aux,inp,i,prec);
    arb_cos(aux,aux,prec);
    arb_neg(aux,aux);
    arb_mul(aux,aux,up->bi+i-1,prec);
    
    //sqrt(k*tanh(k))*k
    arb_set_si(aux2,i);
    arb_tanh(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,i,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_mul_si(aux2,aux2,i,prec);
    
    arb_mul(aux,aux,aux2,prec);
    arb_add(biterm,biterm,aux,prec);
  }
  
  
  arb_add(out,a0kterm,a1kterm,prec);
  arb_add(out,out,biterm,prec);
  
  
  arb_clear(a0kterm);
  arb_clear(a1kterm);
  arb_clear(biterm);
  
  arb_clear(aux); arb_clear(aux2);
  arb_clear(clausen);
  arb_clear(zeta);
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pm_one);
  return 0;
}


// Second primitive of the kernel vanishing at 0 together with its derivative:
// integrating the Fourier series termwise gives
//   IIKp(x) = (1/pi) [ x^2/4 - sum_{m>=1} sqrt(tanh m) (cos(mx) - 1) / m^{5/2} ].
// The sum with weights 1 is Cl_{5/2}(x) - zeta(5/2) (the constant comes from
// the value of the cosine Clausen function at 0); the correction with the
// weights sqrt(tanh m) - 1 is summed up to m = TANH_MAX and its tail is
// bounded by compute_tail_error, using |cos(mx) - 1| <= 2.
void IIKp(arb_t out, const arb_t x, slong prec){
  arb_t argument;
  arb_t error_tail;
  arb_t pi;
  arb_t pm_one;
  arb_t aux,aux2;
  int m;
  
  arb_init(argument);
  arb_init(error_tail);
  arb_init(pi);
  arb_init(pm_one);
  arb_init(aux); arb_init(aux2);
  
  arb_const_pi(pi,prec);
  arb_zero_pm_one(pm_one);
  arb_set_str(argument,"2.5",prec);
  
  // Main term
  arb_clausen_cos(out,argument,x,prec);
  arb_zeta(aux2,argument,prec);
  arb_sub(out,out,aux2,prec);
  
  arb_div_si(aux2,x,2,prec);
  arb_mul(aux2,aux2,aux2,prec);
  arb_sub(out,aux2,out,prec);
  
  // Finite part of the error
  for (m = 1; m <= TANH_MAX; m++){
      
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux2,aux2,1,prec);
      
    arb_mul_si(aux,x,m,prec);
    arb_cos(aux,aux,prec);
    arb_sub_si(aux,aux,1,prec);
    
    arb_mul(aux,aux,aux2,prec);
      
    arb_set_si(aux2,m);
    arb_pow(aux2,aux2,argument,prec);
      
    arb_div(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  
  
  // Tail
  arb_neg(argument,argument);
  compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
  arb_mul_si(error_tail,error_tail,2,prec);
  arb_mul(error_tail,error_tail,pm_one,prec);
  arb_add(out,out,error_tail,prec);
  
  arb_div(out,out,pi,prec);
  
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pi);
  arb_clear(pm_one);
  arb_clear(aux); arb_clear(aux2);
}


void IKp(arb_t out, const arb_t x, slong prec){
  arb_t argument;
  arb_t error_tail;
  arb_t pi;
  arb_t pm_one;
  arb_t aux,aux2;
  int m;
  
  arb_init(argument);
  arb_init(error_tail);
  arb_init(pi);
  arb_init(pm_one);
  arb_init(aux); arb_init(aux2);
  
  arb_const_pi(pi,prec);
  arb_zero_pm_one(pm_one);
  arb_set_str(argument,"1.5",prec);
  
  // Main term
  arb_clausen_sin(out,argument,x,prec);
  arb_div_si(aux2,x,2,prec);
  arb_add(out,out,aux2,prec);
  
  // Finite part of the error
  for (m = 1; m <= TANH_MAX; m++){
      
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux2,aux2,1,prec);
      
    arb_mul_si(aux,x,m,prec);
    arb_sin(aux,aux,prec);
    
    arb_mul(aux,aux,aux2,prec);
      
    arb_set_si(aux2,m);
    arb_pow(aux2,aux2,argument,prec);
      
    arb_div(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  
  // Tail
  arb_neg(argument,argument);
  compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
  arb_mul(error_tail,error_tail,pm_one,prec);
  arb_add(out,out,error_tail,prec);
  
  arb_div(out,out,pi,prec);
  
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pi);
  arb_clear(pm_one);
  arb_clear(aux); arb_clear(aux2);
}

void Kp(arb_t out, const arb_t x, slong prec){
  arb_t argument;
  arb_t error_tail;
  arb_t pi;
  arb_t pm_one;
  arb_t aux,aux2;
  int m;
  
  arb_init(argument);
  arb_init(error_tail);
  arb_init(pi);
  arb_init(pm_one);
  arb_init(aux); arb_init(aux2);
  
  arb_const_pi(pi,prec);
  arb_zero_pm_one(pm_one);
  arb_set_str(argument,"0.5",prec);
  
  // Main term
  arb_clausen_cos(out,argument,x,prec);
  arb_add(out,out,argument,prec);
  
  // Finite part of the error
  for (m = 1; m <= TANH_MAX; m++){
      
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux2,aux2,1,prec);
      
    arb_mul_si(aux,x,m,prec);
    arb_cos(aux,aux,prec);
    
    arb_mul(aux,aux,aux2,prec);
      
    arb_set_si(aux2,m);
    arb_pow(aux2,aux2,argument,prec);
      
    arb_div(aux,aux,aux2,prec);
    arb_add(out,out,aux,prec);
  }
  
  
  // Tail
  arb_neg(argument,argument);
  compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
  arb_mul(error_tail,error_tail,pm_one,prec);
  arb_add(out,out,error_tail,prec);
  
  arb_div(out,out,pi,prec);
  
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pi);
  arb_clear(pm_one);
  arb_clear(aux); arb_clear(aux2);
}


void dKp(arb_t out, const arb_t x, slong prec){
  arb_t argument;
  arb_t error_tail;
  arb_t pi;
  arb_t pm_one;
  arb_t aux,aux2;
  int m;
  
  arb_init(argument);
  arb_init(error_tail);
  arb_init(pi);
  arb_init(pm_one);
  arb_init(aux); arb_init(aux2);
  
  arb_const_pi(pi,prec);
  arb_zero_pm_one(pm_one);
  arb_set_str(argument,"-0.5",prec);
  
  arb_clausen_sin(out,argument,x,prec);
  arb_neg(out,out);
  
  
  // Finite part of the error
  for (m = 1; m <= TANH_MAX; m++){
      
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux2,aux2,1,prec);
      
    arb_mul_si(aux,x,m,prec);
    arb_sin(aux,aux,prec);
    
    arb_mul(aux,aux,aux2,prec);
      
    arb_set_si(aux2,m);
    arb_pow(aux2,aux2,argument,prec);
      
    arb_div(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  
  
  // Tail
  arb_neg(argument,argument);
  compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
  arb_mul(error_tail,error_tail,pm_one,prec);
  arb_add(out,out,error_tail,prec);
  
  arb_div(out,out,pi,prec);
  
  
  
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pi);
  arb_clear(pm_one);
  arb_clear(aux); arb_clear(aux2);
}


void ddKp(arb_t out, const arb_t x, slong prec){
  arb_t argument;
  arb_t error_tail;
  arb_t pi;
  arb_t pm_one;
  arb_t aux,aux2;
  int m;
  
  arb_init(argument);
  arb_init(error_tail);
  arb_init(pi);
  arb_init(pm_one);
  arb_init(aux); arb_init(aux2);
  
  arb_const_pi(pi,prec);
  arb_zero_pm_one(pm_one);
  arb_set_str(argument,"-1.5",prec);
  
  arb_clausen_cos(out,argument,x,prec);
  arb_neg(out,out);
  
  
  // Finite part of the error
  for (m = 1; m <= TANH_MAX; m++){
      
    arb_set_si(aux2,m);
    arb_tanh(aux2,aux2,prec);
    arb_sqrt(aux2,aux2,prec);
    arb_sub_si(aux2,aux2,1,prec);
      
    arb_mul_si(aux,x,m,prec);
    arb_cos(aux,aux,prec);
    
    arb_mul(aux,aux,aux2,prec);
      
    arb_set_si(aux2,m);
    arb_pow(aux2,aux2,argument,prec);
      
    arb_div(aux,aux,aux2,prec);
    arb_sub(out,out,aux,prec);
  }
  
  
  // Tail
  arb_neg(argument,argument);
  compute_tail_error(error_tail,argument,TANH_MAX+1,prec);
  arb_mul(error_tail,error_tail,pm_one,prec);
  arb_add(out,out,error_tail,prec);
  
  arb_div(out,out,pi,prec);
  
  arb_clear(argument);
  arb_clear(error_tail);
  arb_clear(pi);
  arb_clear(pm_one);
  arb_clear(aux); arb_clear(aux2);
  
}


