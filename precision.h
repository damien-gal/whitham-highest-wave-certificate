#ifndef WHITHAM_CAP_PRECISION_H
#define WHITHAM_CAP_PRECISION_H
/* One working precision for every Arb/ACB proof calculation.
   96 bits is insufficient for the retained three-image coefficient test.
   Full-domain validation selects 128 bits; there is no precision escalation. */
#define CAP_PREC 128
#endif
