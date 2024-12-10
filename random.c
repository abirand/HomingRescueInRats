#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "random.h"

void freedist(distribution *d)      // frees dist
{
  free(d->a);
  free(d->p);
  free(d);
}

void sv_rnd(RndState *r, u32b *v)   // save state; also saves for normal, thus v[59]
{
  int i;
  for (i = 0; i < 55; i++)
    v[i] = r->rtab[i];
  v[i] = r->rndx;
  v[i+1] = r->nflg;
  *((double *)(v+i+2)) = r->nv;
}

void rst_rnd(RndState *r, u32b *v)  // restore state; also restores for normal, thus v[59]
{
  int i;

  for (i = 0; i < 55; i++)
    r->rtab[i] = v[i];
  r->rndx = v[i];
  r->nflg = v[i+1];
  r->nv = *((double *)(v+i+2));
}

int nrndm(RndState *r)
{
  int i;

  for (i = 0; i < 24; i++)
    r->rtab[i] -= r->rtab[i+31];
  for (i = 24; i < 55; i++)
    r->rtab[i] -= r->rtab[i-24];
  return 0;
}

void initrand(RndState *r, u32b j)  // Initialize the 32 bit random number generator with seed j
{
  int h, i;
  u32b k;

  for (r->rtab[54] = j |= (k = i = 1); i < 55; i++)
    h = (21*i)%55, r->rtab[--h] = k, k = j - k, j = r->rtab[h];

  while (i--)
    nrndm(r);
  r->rndx = 0;
}

distribution *allocdist(long n)     // allocate a distribution over {0..n-1}
{
  distribution *d;

  d = (distribution *) malloc(sizeof(distribution));
  d->n = n;
  d->a = (u32b *)malloc(d->n * sizeof(long));
  d->p = (double *)malloc(d->n * sizeof(double));
  return d;
}

#define getsmall { while (p[j] >= q) if ((++j) == stop) goto end; t = j++; }
#define getlarge   while (p[k] <  q) if ((++k) == stop) goto cleanup;

distribution *initdist(distribution *d, double s) // Initialize the distribution d
{
  /*
  Note: d->p must have d->n elements which sum to s on entry to initdist.
 The elements of d->p and d->a are overwritten by the initialization process.
 */

  u32b j, k, t, stop, *a;
  double q, *p;

  stop = d->n, q = s/stop, j = k = 0;

  d->m1 = stop/TWO_LONGSZ;
  d->m2 = s/(stop * TWO_LONGSZ);

  a = d->a;
  p = d->p;

  getsmall;
  getlarge;

  loop:

  a[t] = k;
  p[k] += p[t] - q;

  if (p[k] >= q) {
    if (j == stop)
      goto end;
    getsmall;
    goto loop;
  }
  t = k++;
  if (k == stop)
    goto cleanup;
  if (j < k)getsmall;
    getlarge;
  goto loop;

  cleanup:

  a[t] = t;
  while( j < stop )
  {
    a[j] = j;
    j++;
  }

  end: return d;
}

u32b drand(RndState *r, distribution *d)  // Return element from {0..d->n-1} according to d
{
  int j;
  u32b r1 = rndm(r);

  if ((rndm(r)+(r1&0xffff)/65536.0)*(d->m2) < d->p[j=r1*(d->m1)])
    return j;
  return d->a[j];
}

double normal(RndState *rn, double m, double s) // normaly distributed; mean m standard deviation s
{
  // nflg nd nv are globals (so state can easily be saved/restored)
  double u, r;

  if( !rn->nflg )
  {
    do
      u = 2*U01(rn) - 1., rn->nv = 2*U01(rn) - 1.;
    while ((r = u*u + rn->nv*rn->nv) >= 1 || (r == 0));
    r = sqrt( -2.*log(r)/r);
    u *= r;
    rn->nv *= r;
    rn->nflg = 1;
    return s*u + m;
  }
  else
  {
    rn->nflg = 0;
    return s*rn->nv + m;
  }
}

u32b Poisson(RndState *r, double mu)
{
  double p = exp( -mu ), q = 1;
  int    n = -1;

  do q *= U01(r), n++;
  while( (q > p) || (q == p) );

  return n;
}
