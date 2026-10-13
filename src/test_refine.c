/* Independent exhaustive checks of the discrete convex minimization against
 * its literal two-variable optimization on small integer domains. */
#include "refine.c"
#include <stdio.h>
static unsigned state=1234567;
static unsigned rnd(void){state=1664525*state+1013904223;return state;}
int main(void){
 for(int trial=0;trial<300;trial++){
  V a[8],b[8],m[8],v[8];v[0]=0;
  for(int j=0;j<8;j++){a[j]=rnd()%10000;b[j]=rnd()%10000;m[j]=rnd()%10000;if(j)v[j]=v[j-1]+rnd()%8;}
  Row r={0};r.difference=1;r.a=make(a,v,8);r.b=make(b,v,8);r.m=make(m,v,8);r.pu=minp(&r,0,48);r.tu=mint(&r,0,64);
  for(V t=0;t<=64;t++){
   I best=0;int first=1;
   for(V p=0;p<=48&&p<=t;p++){
    I val=square(p);for(int j=0;j<8;j++)val-=(I)b[j]*(v[j]<p?v[j]:p);
    for(int j=0;j<8;j++)val+=(I)a[j]*(v[j]>t?v[j]:t)-(I)m[j]*(v[j]<t?v[j]:t);
    if(first||val<best){first=0;best=val;}
   }
   if(inner(&r,t)!=best){fprintf(stderr,"FAIL inner minimization\n");return 1;}
  }
  for(V q=0;q<=64;q++){
   I best=inner(&r,0);for(V t=1;t<=q;t++)if(inner(&r,t)<best)best=inner(&r,t);
   if(inner(&r,q<r.tu?q:r.tu)!=best){fprintf(stderr,"FAIL clipped minimizer\n");return 1;}
  }
  clear(&r.a);clear(&r.b);clear(&r.m);
 }
 puts("PASS 300 exhaustive discrete-convex problems, all clipped domains");return 0;
}
