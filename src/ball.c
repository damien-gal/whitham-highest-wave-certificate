#include "ball.h"
#include <math.h>
#define BLOCK 4096
#define MAX_BLOCKS 4096
static arb_ptr blocks[MAX_BLOCKS];static size_t used;
void require(int ok,const char *why){if(!ok){fprintf(stderr,"UNRESOLVED: %s\n",why);exit(2);}}
size_t mark(void){return used;}
arb_ptr tmp(void){size_t k=used/BLOCK;require(k<MAX_BLOCKS,"scratch arena capacity");if(!blocks[k])blocks[k]=_arb_vec_init(BLOCK);return blocks[k]+used++%BLOCK;}
void reset(size_t m){require(m<=used,"scratch scope order");while(used>m){--used;arb_zero(blocks[used/BLOCK]+used%BLOCK);}}
arb_ptr integer(slong n){arb_ptr z=tmp();arb_set_si(z,n);return z;}
arb_ptr decimal(const char *s){arb_ptr z=tmp();require(!arb_set_str(z,s,CAP_PREC),"decimal input");return z;}
arb_ptr rational(slong n,ulong d){arb_ptr z=tmp();fmpq_t q;fmpq_init(q);require(d>0,"rational denominator");fmpq_set_si(q,n,d);arb_set_fmpq(z,q,CAP_PREC);fmpq_clear(q);return z;}
#define BIN(name,op) arb_ptr name(B a,B b){arb_ptr z=tmp();op(z,a,b,CAP_PREC);return z;}
BIN(add,arb_add) BIN(sub,arb_sub) BIN(mul,arb_mul) BIN(divide,arb_div)
BIN(minimum_b,arb_min) BIN(maximum_b,arb_max)
#define UN(name,op) arb_ptr name(B a){arb_ptr z=tmp();op(z,a,CAP_PREC);return z;}
UN(root,arb_sqrt) UN(sq,arb_sqr) UN(sinb,arb_sin) UN(cosb,arb_cos) UN(expb,arb_exp) UN(zetab,arb_zeta) UN(asinb,arb_asin) UN(asinhb,arb_asinh) UN(acoshb,arb_acosh)
arb_ptr neg(B a){arb_ptr z=tmp();arb_neg(z,a);return z;}
arb_ptr ab(B a){arb_ptr z=tmp();arb_abs(z,a);return z;}
arb_ptr nonneg(B a){arb_ptr z=tmp();arb_nonnegative_part(z,a);return z;}
arb_ptr lower(B a){arb_ptr z=tmp();arb_get_lbound_arf(arb_midref(z),a,CAP_PREC);mag_zero(arb_radref(z));return z;}
arb_ptr upper(B a){arb_ptr z=tmp();arb_get_ubound_arf(arb_midref(z),a,CAP_PREC);mag_zero(arb_radref(z));return z;}
arb_ptr absupper(B a){return upper(ab(a));}
arb_ptr hull(B a,B b){arb_ptr z=tmp();arb_union(z,a,b,CAP_PREC);return z;}
arb_ptr inflate(B a,B e){arb_ptr z=tmp();arb_set(z,a);arb_add_error(z,e);return z;}
arb_ptr power(B a,ulong n){arb_ptr z=tmp();arb_pow_ui(z,a,n,CAP_PREC);return z;}
arb_ptr times(B a,slong n){arb_ptr z=tmp();arb_mul_si(z,a,n,CAP_PREC);return z;}
arb_ptr over(B a,ulong n){arb_ptr z=tmp();arb_div_ui(z,a,n,CAP_PREC);return z;}
arb_ptr dyadic(int64_t n,int bits){arb_ptr z=integer(n);arb_mul_2exp_si(z,z,-bits);return z;}
arb_ptr exact_double(double x){require(isfinite(x),"finite binary64 input");arb_ptr z=tmp();arb_set_d(z,x);return z;}
int positive(B a){return arb_is_positive(a);}int negative(B a){return arb_is_negative(a);}int lt(B a,B b){return arb_lt(a,b);}int le(B a,B b){return arb_le(a,b);}
int64_t lattice(B a,int bits,int up){size_t m=mark();require(arb_is_finite(a),"finite lattice conversion");B v=up?upper(a):lower(a);arb_ptr z=tmp();arb_mul_2exp_si(z,v,bits);fmpz_t n;fmpz_init(n);arf_get_fmpz(n,arb_midref(z),up?ARF_RND_CEIL:ARF_RND_FLOOR);require(fmpz_fits_si(n),"lattice fits signed word");int64_t out=fmpz_get_si(n);fmpz_clear(n);reset(m);return out;}
void show(const char *name,B a){require(arb_is_finite(a),name);printf("%s = ",name);arb_printn(a,25,0);puts("");char *s=arb_dump_str(a);printf("ARB_DUMP %s: %s\n",name,s);flint_free(s);fflush(stdout);}
