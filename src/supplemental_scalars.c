#include "constant_functions.h"
int main(void) {
 const slong p=CAP_PREC; u_params_struct up; populate_parameters(&up,p);
 arb_t e,v,t,cb,a,b,pi,z;arb_init(e);arb_init(v);arb_init(t);arb_init(cb);arb_init(a);arb_init(b);arb_init(pi);arb_init(z);
 arb_set_str(e,"0.1",p);int ok=1;
 for(int M=5;M<=7;M+=2){
 printf("M=%d\n",M);
 compute_c_u0_hat_minus(v,e,M,&up,p); arb_set_str(t,".217",p);ok&=arb_gt(v,t);printf("hat lower ");arb_printn(v,25,0);puts("");
 compute_c_u0_hat_plus(v,e,M,&up,p); arb_set_str(t,".24",p);ok&=arb_lt(v,t);if(M==5){arb_set_str(t,".2395495",p);ok&=arb_lt(v,t);}printf("hat upper ");arb_printn(v,25,0);puts("");
 compute_c_p0(v,e,M,&up,p); arb_set_str(t,".142",p);ok&=arb_lt(v,t);if(M==5){arb_set_str(t,".1418003",p);ok&=arb_lt(v,t);}printf("cp0 upper ");arb_printn(v,25,0);puts("");
 }
 /* Printed local T0 criterion with rounded manuscript c_p0=.142. */
 arb_set_str(e,".05",p);arb_set_str(cb,".997362",p);arb_set_str(v,".142",p);arb_mul(v,v,cb,p);arb_set_str(z,".5",p);arb_sub(z,up.p0,z,p);arb_pow(z,e,z,p);arb_mul(v,v,z,p);
 compute_cT0_1(a,e,&up,p);compute_cT0_2(b,e,&up,p);arb_add(a,a,b,p);arb_sqr(z,e,p);arb_mul(a,a,z,p);arb_add(v,v,a,p);
 arb_const_pi(pi,p);arb_set_ui(a,2);arb_sqrt(a,a,p);arb_sub_ui(a,a,2,p);arb_div_ui(a,a,2,p);arb_div(a,a,pi,p);arb_set_str(b,".5",p);arb_zeta(b,b,p);arb_mul(a,a,b,p);arb_sub(v,v,a,p);
 ok&=arb_is_negative(v);arb_set_str(t,"-.034",p);ok&=arb_lt(v,t); /* Appendix B prints < -0.034 */
 printf("printed local T0 residual ");arb_printn(v,25,0);puts("");puts(ok?"PASS supplemental scalar comparisons":"UNRESOLVED");return !ok;
}
