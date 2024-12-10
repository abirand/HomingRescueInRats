#ifndef   RANDOM__
#define   RANDOM__

#include "types.h"

#define   TWO_32     (4294967296.0)
#define   TWO_LONGSZ TWO_32

typedef struct RndState
{
  u32b    rtab[55];
  int     rndx,nflg;
  double  nv;
} RndState;

#define rndm(r)   ((++(r)->rndx>54)?(r)->rtab[(r)->rndx=nrndm(r)]:(r)->rtab[(r)->rndx]) /* random 32-bit generator  */
#define U01(r)    (rndm(r)/TWO_LONGSZ)                            /* random in interval [0,1) */
#define U(r,x)    (U01(r)*(x))                                    /* random in interval [0,x) */
#define rnd(r,n)  ((u32b)U(r,n))                                  /* random from set {0..n-1} */

// for generating random number i with probability p[i]
typedef struct dist
{
  double *p;
  u32b   *a;
  u32b    n;
  double  m1;
  double  m2;
} distribution;

#ifdef __cplusplus
extern "C"
{
#endif

void          sv_rnd(RndState *r, u32b *v);             // save state; also saves for normal, thus v[59]
void          rst_rnd(RndState *r, u32b *v);            // restore state; also restores for normal, thus v[59]
int           nrndm(RndState *r);
void          initrand(RndState *r, u32b j);            // Initialize the 32 bit random number generator with seed j
distribution *allocdist(long n);                        // allocate a distribution over {0..n-1}
distribution *initdist(distribution *d, double s);      // Initialize the distribution d
void          freedist(distribution *d);                // frees dist
u32b          drand(RndState *r, distribution *d);      // Return element from {0..d->n-1} according to d
double        normal(RndState *r, double m, double s);

// Returns choice from Poisson distribution with parameter (mean) mu
u32b          Poisson(RndState *r, double mu);

#ifdef __cplusplus
}
#endif

#endif 
/* RANDOM__ */
