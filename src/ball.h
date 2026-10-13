#ifndef WHITHAM_BALL_H
#define WHITHAM_BALL_H
/* Single-thread scratch arena for readable Arb formulas. Every operation
   allocates a fresh initialized ball; reset(mark) releases only later balls.
   Persistent values must be copied with arb_set before resetting a scope.
   No expression aliases or floating-point approximations are used. */
#include <arb.h>
#include <fmpq.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef CAP_PREC
#define CAP_PREC 128
#endif
typedef arb_srcptr B;
void require(int ok,const char *why);
size_t mark(void);void reset(size_t m);arb_ptr tmp(void);
arb_ptr integer(slong n);arb_ptr decimal(const char *s);arb_ptr rational(slong n,ulong d);
arb_ptr add(B a,B b);arb_ptr sub(B a,B b);arb_ptr mul(B a,B b);arb_ptr divide(B a,B b);
arb_ptr neg(B a);arb_ptr ab(B a);arb_ptr sq(B a);arb_ptr root(B a);arb_ptr power(B a,ulong n);
arb_ptr sinb(B a);arb_ptr cosb(B a);arb_ptr expb(B a);arb_ptr zetab(B a);arb_ptr asinb(B a);arb_ptr asinhb(B a);arb_ptr acoshb(B a);
arb_ptr lower(B a);arb_ptr upper(B a);arb_ptr absupper(B a);arb_ptr hull(B a,B b);arb_ptr minimum_b(B a,B b);arb_ptr maximum_b(B a,B b);arb_ptr nonneg(B a);arb_ptr inflate(B a,B e);
arb_ptr times(B a,slong n);arb_ptr over(B a,ulong n);arb_ptr dyadic(int64_t n,int bits);arb_ptr exact_double(double x);
int positive(B a);int negative(B a);int lt(B a,B b);int le(B a,B b);
int64_t lattice(B a,int bits,int up);
void show(const char *name,B a);
#endif
