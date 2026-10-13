/* Direct C/Arb port of the normalized kernel-zero certificate. Closed cells
   are covered recursively: both children must pass. All targets are unchanged. */
#include "kernel.h"
#include <string.h>
static Kernel k;static arb_ptr coeff[2],gap;static long leaves,nodes;static int reject_test;
static int cell(B a,B b,int depth,int cusp){
 size_t m=mark();nodes++;arb_ptr v[2],d[2],c[2];for(int s=0;s<2;s++){v[s]=tmp();d[s]=tmp();c[s]=tmp();poly_box(v[s],d[s],c[s],coeff[s],13,a,b);}
 arb_ptr gv=tmp(),gd=tmp(),gc=tmp();poly_box(gv,gd,gc,gap,13,a,b);
 int ok=positive(v[0])&&lt(v[1],integer(1))&&arb_is_nonnegative(gv)&&le(gv,decimal(".000013"));B w[2]={NULL,NULL};
 if(ok){
  if(cusp){ok=le(decimal(".652894"),v[0])&&le(v[0],decimal(".652897"))&&le(decimal(".652904"),v[1])&&le(v[1],decimal(".652907"));for(int s=0;s<2;s++)w[s]=small_curve(&k,v[s],b);}
  else {B mid=over(add(a,b),2),X=hull(a,b),h=inflate(integer(0),upper(over(sub(b,a),2)));for(int s=0;s<2;s++)w[s]=add(curve(&k,mid,c[s]),mul(curve_derivative(&k,X,v[s],d[s]),h));}
  B bound=maximum_b(mul(absupper(w[0]),root(divide(v[1],v[0]))),absupper(w[1]));
  ok=ok&&negative(w[0])&&positive(w[1])&&le(bound,decimal(reject_test?"0":".000012"));
 }
 if(ok){leaves++;printf("SIGN_LEAF %ld depth=%d cusp=%d ",leaves,depth,cusp);arb_printn(a,18,0);putchar(' ');arb_printn(b,18,0);puts("");reset(m);return 1;}
 if(depth>=15||reject_test){reset(m);return 0;}
 B mid=over(add(a,b),2);ok=cell(a,mid,depth+1,cusp)&&cell(mid,b,depth+1,cusp);reset(m);return ok;
}
int main(int argc,char **argv){reject_test=argc>1&&!strcmp(argv[1],"reject");kernel_init(&k,80,40);coeff[0]=_arb_vec_init(13);coeff[1]=_arb_vec_init(13);gap=_arb_vec_init(13);root_coefficients(coeff[0],coeff[1]);for(int i=0;i<13;i++)arb_sub(gap+i,coeff[1]+i,coeff[0]+i,CAP_PREC);
 require(lt(k.pi,rational(22,7)),"rational cover beyond pi");
 for(int i=0;i<64;i++){size_t m=mark();require(cell(rational(i,6400),rational(i+1,6400),0,1),"cusp sign cover");reset(m);}
 for(int i=0;i<384;i++){size_t m=mark();B eps=rational(1,100),span=sub(rational(22,7),eps);require(cell(add(eps,mul(span,rational(i,384))),add(eps,mul(span,rational(i+1,384))),0,0),"compact sign cover");reset(m);}
 printf("PASS kernel-sign full [0,22/7] precision=%d leaves=%ld nodes=%ld\n",CAP_PREC,leaves,nodes);return 0;
}
