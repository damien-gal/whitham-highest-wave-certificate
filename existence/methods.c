/* Blah blah blah. Javi */

#include "basic_functions.h"

/* TANH_MAX comes from basic_functions.h */
#define QSIZE 1024


// Checks if the function "func" is smaller than the bound "bound".
// We subdivide until the minimum_width "min_width".
int compute_bound_C0(const arb_t min_width, const arb_t bound, whitham_calc_func func, void * param, const arb_t inp, slong prec){
  arb_t curr;
  arb_t rcurr;
  arb_t left,right,mid;
  arb_t fl,fr;
  arf_t arfleft,arfright;
  
  int bl, br;
    
  arb_init(curr);
  arb_init(rcurr);
  arb_init(left); arb_init(right); arb_init(mid);
  arb_init(fl); arb_init(fr);
  arf_init(arfleft); arf_init(arfright);
  
  func(curr,inp,param,prec);
  arb_get_rad_arb(rcurr,inp);
  printf("%s  %s\n",arb_get_str(inp,20,1),arb_get_str(curr,20,1));
  if (arb_lt(curr,bound)){
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr);
    arf_clear(arfleft); arf_clear(arfright);
    return 1;
  }
  if (arb_gt(curr,bound)){
    printf("Violated bound at x = ");
    arb_printd(inp,prec);
    printf("\n");
    printf("f(x) = ");
    arb_printd(curr,prec);
    printf("\n");
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  if (arb_lt(rcurr,min_width)){
    printf("Attained minimum_width\n");
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  arb_get_mid_arb(mid,inp);
  arb_get_interval_arf(arfleft, arfright, inp, prec);
  arb_set_arf(left,arfleft);
  arb_set_arf(right,arfright);
  arb_union(left,left,mid,prec);
  arb_union(right,mid,right,prec);
    
  bl = compute_bound_C0(min_width,bound,func,param,left,prec);
  if (!bl){
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  br = compute_bound_C0(min_width,bound,func,param,right,prec);
  
  arb_clear(curr);
  arb_clear(rcurr);
  arb_clear(left); arb_clear(right); arb_clear(mid);
  arb_clear(fl); arb_clear(fr);
  arf_clear(arfleft); arf_clear(arfright);
  return bl & br;
}



// Checks if the function "func" is smaller than the bound "bound".
// We subdivide until the minimum_width "min_width".
// We use C^1 information
int compute_bound_C1(const arb_t min_width, const arb_t bound, whitham_calc_func func, whitham_calc_func deriv, void * param, const arb_t inp, slong prec){
  arb_t curr, currc0;
  arb_t dcurr;
  arb_t rcurr;
  arb_t left,right,mid;
  arb_t fl,fr,fmid;
  arf_t arfleft,arfright;
  
  int bl, br;
    
  arb_init(curr); arb_init(currc0);
  arb_init(dcurr);
  arb_init(rcurr);
  arb_init(left); arb_init(right); arb_init(mid);
  arb_init(fl); arb_init(fr); arb_init(fmid);
  arf_init(arfleft); arf_init(arfright);
  
  arb_get_mid_arb(mid,inp);
  
  func(fmid,mid,param,prec);
  deriv(dcurr,inp,param,prec);
  
  arb_sub(curr,inp,mid,prec);
  arb_mul(curr,curr,dcurr,prec);
  arb_add(curr,curr,fmid,prec);
    
  arb_get_rad_arb(rcurr,inp);
  
  printf("%s  %s\n",arb_get_str(inp,20,1),arb_get_str(curr,20,1));
  if (arb_lt(curr,bound)){
    arb_clear(curr); arb_clear(currc0);
    arb_clear(dcurr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr); arb_clear(fmid);
    arf_clear(arfleft); arf_clear(arfright);
    return 1;
  }
  if (arb_gt(curr,bound)){
    printf("Violated bound at x = ");
    arb_printd(inp,prec);
    printf("\n");
    printf("f(x) = ");
    arb_printd(curr,prec);
    printf("\n");
    arb_clear(curr); arb_clear(currc0);
    arb_clear(dcurr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr); arb_clear(fmid);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  if (arb_lt(rcurr,min_width)){
    printf("Attained minimum_width\n");
    arb_clear(curr); arb_clear(currc0);
    arb_clear(dcurr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr); arb_clear(fmid);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  arb_get_interval_arf(arfleft, arfright, inp, prec);
  arb_set_arf(left,arfleft);
  arb_set_arf(right,arfright);
  arb_union(left,left,mid,prec);
  arb_union(right,mid,right,prec);
    
  bl = compute_bound_C1(min_width,bound,func,deriv,param,left,prec);
  if (!bl){
    arb_clear(curr); arb_clear(currc0);
    arb_clear(dcurr);
    arb_clear(rcurr);
    arb_clear(left); arb_clear(right); arb_clear(mid);
    arb_clear(fl); arb_clear(fr); arb_clear(fmid);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  br = compute_bound_C1(min_width,bound,func,deriv,param,right,prec);
  
  arb_clear(curr); arb_clear(currc0);
  arb_clear(dcurr);
  arb_clear(rcurr);
  arb_clear(left); arb_clear(right); arb_clear(mid);
  arb_clear(fl); arb_clear(fr); arb_clear(fmid);
  arf_clear(arfleft); arf_clear(arfright);
  return bl & br;
}

// TODO(Javi): untested
// TODO(Javi): make some control to return 0 is something fails
int MidPointQuadratureAbs(arb_t out, const arb_t inp_left, const arb_t inp_right, void* param, 
                       whitham_calc_func Integrand_Function, whitham_calc_func Integrand_Derivative, slong prec){
  arb_t uinp;
  arb_t dx;
  
  arb_init(uinp);
  arb_init(dx);
  
  arb_sub(dx,inp_right,inp_left,prec);
  arb_union(uinp,inp_left,inp_right,prec);
  
  arb_t mid; arb_t fmid; arb_t d2fmid; arb_t aux;
  arb_init(mid); arb_init(fmid); arb_init(d2fmid); arb_init(aux);
  
  arb_add(mid,inp_left,inp_right,prec); arb_div_si(mid,mid,2,prec);
  
  Integrand_Function(out,uinp,param,prec);
  arb_mul(out,out,dx,prec);  
  
  if (arb_is_nonpositive(out) || arb_is_nonnegative(out)){ // Sign condition, use midpoint quadrature
    // TODO(Javi): Perhaps implement error control below the next 2 lines?
    Integrand_Function(fmid,mid,param,prec);
    Integrand_Derivative(d2fmid,uinp,param,prec);
    arb_div_si(aux,d2fmid,24,prec);
    arb_mul(aux,aux,dx,prec);
    arb_mul(aux,aux,dx,prec);
    arb_add(aux,aux,fmid,prec);
    arb_mul(aux,aux,dx,prec);
    arb_intersection(out,out,aux,prec); // Take the intersection between C^0 and C^2 enclosures
  }
  else { // No sign condition, use C^0 enclosure
  }
  
  arb_abs(out,out);
  
  
  arb_clear(uinp);
  arb_clear(dx);
  arb_clear(mid); arb_clear(fmid); arb_clear(d2fmid); arb_clear(aux);
  return 1;
}

int integrate(arb_t out, const arb_t inp_left, const arb_t inp_right, void* param, 
              whitham_calc_func Integrand_Function, whitham_calc_func Integrand_Derivative, slong prec,
              const arb_t max_width_abs, const arb_t max_width_rel, int error_mask){
  arb_ptr queue_regions_left,queue_regions_right;
  arb_t curr_left,curr_right;
  arb_t fcurr;
  arb_t rcurr,rcurr_rel_th;
  arb_t mid;
  
  int start;
  int end;
  
  queue_regions_left = _arb_vec_init(QSIZE);
  queue_regions_right = _arb_vec_init(QSIZE);
  arb_init(curr_left); arb_init(curr_right);
  arb_init(fcurr);
  arb_init(rcurr); arb_init(rcurr_rel_th);
  arb_init(mid);
  
  arb_set(queue_regions_left,inp_left);
  arb_set(queue_regions_right,inp_right);
  
  arb_zero(out);
  
  start = 0;
  end = 1;
  
  
  while (start != end){
    arb_set(curr_left,queue_regions_left+start);
    arb_set(curr_right,queue_regions_right+start);
    start++; start%= QSIZE;
    if (!MidPointQuadratureAbs(fcurr, curr_left, curr_right, param, Integrand_Function, Integrand_Derivative, prec)) return 0;
    arb_get_rad_arb(rcurr,fcurr);
    
    // rcurr_rel_th = dx*max_width_rel.
    arb_sub(rcurr_rel_th,curr_right,curr_left,prec);
    arb_mul(rcurr_rel_th,rcurr_rel_th,max_width_rel,prec);
    
    int pass_error = 1;
    if ((error_mask & 1) && !arb_lt(rcurr,max_width_abs)) pass_error = 0;
    if (((error_mask >> 1) & 1) && !arb_lt(rcurr,rcurr_rel_th)) pass_error = 0;
    
    if (pass_error){
      arb_add(out,out,fcurr,prec);
    }
    else { //split
      
      arb_add(mid,curr_left,curr_right,prec);
      arb_div_si(mid,mid,2,prec);
            
      arb_set(queue_regions_left+end,curr_left);
      arb_set(queue_regions_right+end,mid);
      end++; end%=QSIZE;
      if (end == start){
        printf("Exceeded Size of Queue: %d\n",QSIZE);
        arb_clear(queue_regions_left);
        arb_clear(queue_regions_right);
        arb_clear(curr_left); arb_clear(curr_right);
        arb_clear(fcurr);
        arb_clear(rcurr); arb_clear(rcurr_rel_th);
        arb_clear(mid);
        return 0;
      }
      arb_set(queue_regions_left+end,mid);
      arb_set(queue_regions_right+end,curr_right);
      end++; end%=QSIZE;
      if (end == start){
        printf("Exceeded Size of Queue: %d\n",QSIZE);
        arb_clear(queue_regions_left);
        arb_clear(queue_regions_right);
        arb_clear(curr_left); arb_clear(curr_right);
        arb_clear(fcurr);
        arb_clear(rcurr); arb_clear(rcurr_rel_th);
        arb_clear(mid);
        return 0;
      }
    }
  }
  
  
  arb_clear(queue_regions_left);
  arb_clear(queue_regions_right);
  arb_clear(curr_left); arb_clear(curr_right);
  arb_clear(fcurr);
  arb_clear(rcurr); arb_clear(rcurr_rel_th);
  arb_clear(mid);
  printf("Everything okay!\n");
  return 1;
}


// TODO(Javi): OJO. inp_left, inp_right son valores de w, no de x!!!
int bound_operator(const arb_t inp_left, const arb_t inp_right, void* param, 
                   whitham_calc_func Integrand_Function, whitham_calc_func Integrand_Derivative,
                   explicit_integral BigY_Function, explicit_singularity SingularityZero_Function,
                   explicit_singularity SingularityOne_Function, const arb_t delta0, const arb_t delta1,
                   slong prec, const arb_t min_width_x, const arb_t max_width_abs, const arb_t max_width_rel, 
                   int error_mask, const arb_t bound){
  Integration_params_struct* Ips = param;
  arb_t curr;
  arb_t rcurr;
  arb_t xleft,xright,xmid;
  arb_t fl,fr;
  arb_t big_y;
  arb_t singzero, singone;
  arf_t arfleft,arfright;
  
  int bl, br;
  int outcome_integration;
    
  arb_init(curr);
  arb_init(rcurr);
  arb_init(xleft); arb_init(xright); arb_init(xmid);
  arb_init(fl); arb_init(fr);
  arb_init(big_y);
  arb_init(singzero); arb_init(singone);
  arf_init(arfleft); arf_init(arfright);
   
  
  // Main term: y < x
  outcome_integration = integrate(curr, inp_left, inp_right, param, Integrand_Function, Integrand_Derivative, prec, max_width_abs,
                                  max_width_rel, error_mask);

  // y > x
  if (outcome_integration){
    outcome_integration = outcome_integration & BigY_Function(big_y,param,prec);
  }
  // Singularities
  if (outcome_integration){
    outcome_integration = outcome_integration & SingularityZero_Function(singzero,delta0,param,prec);
    outcome_integration = outcome_integration & SingularityOne_Function(singone,delta1,param,prec);
  }
  
  arb_add(curr,curr,big_y,prec);
  arb_add(curr,curr,singzero,prec);
  arb_add(curr,curr,singone,prec);
  
  arb_get_rad_arb(rcurr,Ips->x);
  
  printf("%s  %s %s\n",arb_get_str(Ips->x,20,1),arb_get_str(curr,20,1),arb_get_str(big_y,20,1));
  if (outcome_integration && arb_lt(curr,bound)){
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(xleft); arb_clear(xright); arb_clear(xmid);
    arb_clear(fl); arb_clear(fr);
    arb_clear(big_y);
    arf_clear(arfleft); arf_clear(arfright);
    return 1;
  }
  if (outcome_integration && arb_gt(curr,bound)){
    printf("Violated bound at x = ");
    arb_printd(Ips->x,prec);
    printf("\n");
    printf("f(x) = ");
    arb_printd(curr,prec);
    printf("\n");
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(xleft); arb_clear(xright); arb_clear(xmid);
    arb_clear(fl); arb_clear(fr);
    arb_clear(big_y);
    arb_clear(singzero); arb_clear(singone);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  if (arb_lt(rcurr,min_width_x)){
    printf("rcurr = %s min width = %s\n",arb_get_str(rcurr,20,1),arb_get_str(min_width_x,20,1));
    printf("Attained minimum_width\n");
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(xleft); arb_clear(xright); arb_clear(xmid);
    arb_clear(fl); arb_clear(fr);
    arb_clear(big_y);
    arb_clear(singzero); arb_clear(singone);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  arb_get_mid_arb(xmid,Ips->x);
  arb_get_interval_arf(arfleft, arfright, Ips->x, prec);
  arb_set_arf(xleft,arfleft);
  arb_set_arf(xright,arfright);
  arb_union(xleft,xleft,xmid,prec);
  arb_union(xright,xmid,xright,prec);
    
  populate_Integration_parameters(Ips, xleft, Ips->order, prec);
  
  bl = bound_operator(inp_left, inp_right, Ips, Integrand_Function, Integrand_Derivative, BigY_Function,
                      SingularityZero_Function, SingularityOne_Function, delta0, delta1, prec, min_width_x, max_width_abs, max_width_rel, error_mask, bound);
                      
  if (!bl){
    arb_clear(curr);
    arb_clear(rcurr);
    arb_clear(xleft); arb_clear(xright); arb_clear(xmid);
    arb_clear(fl); arb_clear(fr);
    arb_clear(big_y);
    arb_clear(singzero); arb_clear(singone);
    arf_clear(arfleft); arf_clear(arfright);
    return 0;
  }
  populate_Integration_parameters(Ips, xright, Ips->order, prec);
  
  br = bound_operator(inp_left, inp_right, Ips, Integrand_Function, Integrand_Derivative, BigY_Function,
                      SingularityZero_Function, SingularityOne_Function, delta0, delta1, prec, min_width_x, max_width_abs, max_width_rel, error_mask, bound);
    
  arb_clear(curr);
  arb_clear(rcurr);
  arb_clear(xleft); arb_clear(xright); arb_clear(xmid);
  arb_clear(fl); arb_clear(fr);
  arb_clear(big_y);
  arb_clear(singzero); arb_clear(singone);
  arf_clear(arfleft); arf_clear(arfright);
  return bl & br;
}

