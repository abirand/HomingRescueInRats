#include <math.h>

// Useful "utility" functions
// Similar to fprintf(stderr...); exit(-1);
static void Error(const char *fmt, ...)
{
  va_list ap;

  if (!fmt) return;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}


// This is Michael's version of c's exp() function.
// I don't know much about it except it was used in the sim project,
// and is supposed to be a good bit faster than c's version.
static double Exp(double x)
{
  int i,j;  double p;  static double T[512];  static int t = 0;

  if( !t ) {
    // initialize Exp array
    for (i = 512; i--; T[i] = exp((double)i)); 
    t = 1;
  }

  if (x > 0){
    if ((i = x) > 511)
      return 2.2844135865397566e+222;
    for (p = 1. + (x-i)/256., j = 8; j--; p *= p);
    return p*T[i];
  }
  if (x < 0){
    if ((i = -x) > 511)
      return 4.377491037053051e-223;
    for (p = 1. - (x+i)/256., j = 8; j--; p *= p);
    return 1./(p*T[i]);
  }
  return 1.;
}

// Re-mapping macro
#if _FAST_EXP
#define exp(x) Exp(x)
#endif


// Functions dealing directly with maintenance of Individuals
// Returns a pointer to an available individual structure.  NULL for none.
static pIndividual GetCurrentIndividual(pThreadData td)
{
  // Record where this individual was allocated from
  td->ICurrent[td->nICurrent]->id = td->nICurrent;
  return td->ICurrent[td->nICurrent++];
}

// Links an individual's x, y, and z characteristics to it's data (d) member
static void LinkIndividual(pIndividual i)
{

  i->x0 = i->d;
  i->x1 = i->d+LENGTH_E;
  i->m0 = i->d+(LENGTH_E<<1);
  i->m1 = i->d+(LENGTH_E<<1)+LENGTH_M;
  i->k0 = i->d+(LENGTH_E<<1)+(LENGTH_M<<1);
  i->k1 = i->d+(LENGTH_E<<1)+(LENGTH_M<<1)+LENGTH_K;
  i->f0 = i->d+(LENGTH_E<<1)+(LENGTH_M<<1)+(LENGTH_K<<1);
  i->f1 = i->d+(LENGTH_E<<1)+(LENGTH_M<<1)+(LENGTH_K<<1)+LENGTH_F;
  i->z0 = i->d+(LENGTH_E<<1)+(LENGTH_M<<1)+(LENGTH_K<<1)+(LENGTH_F<<1);
  i->z1 = i->d+(LENGTH_E<<1)+(LENGTH_M<<1)+(LENGTH_K<<1)+(LENGTH_F<<1)+LENGTH_N;
}


/*
  Sets the locus l in allele array c to the value v
  Individuals have 4 allele arrays:
  ->x0  // Characteristic alleles 1
  ->x1  // Characteristic alleles 2
  
  w stands for 'which'.  This identifies either:
  ecological traits: w == 0
  mating     traits: w == 1
  marker     traits: w == 2
  fpref      traits: w == 3
  neutral    traits: w == 4

  To set a specific locus to a specific value, use SetLocus() like so:
  SetLocus(i->x0, 3, 0);
  This will set the third locus in the first characteristic allele array to 0
*/
static void SetLocus(u64b *c, u64b l, u64b v, int w)
{
  // Clear out the individual's current locus and or-in the new locus
  switch(w) {
  case 0:
    c[l/LPW_E] &= ~(MASK_E<<(_Ae*(l%LPW_E)));
    c[l/LPW_E] |= (v<<(_Ae*(l%LPW_E))) & (MASK_E<<(_Ae*(l%LPW_E)));
    break;
  case 1:
    c[l/LPW_M] &= ~(MASK_M<<(_Am*(l%LPW_M)));
    c[l/LPW_M] |= (v<<(_Am*(l%LPW_M))) & (MASK_M<<(_Am*(l%LPW_M)));
    break;
  case 2:
    c[l/LPW_K] &= ~(MASK_K<<(_Ak*(l%LPW_K)));
    c[l/LPW_K] |= (v<<(_Ak*(l%LPW_K))) & (MASK_K<<(_Ak*(l%LPW_K)));
    break;
  case 3:
    c[l/LPW_F] &= ~(MASK_F<<(_Af*(l%LPW_F)));
    c[l/LPW_F] |= (v<<(_Af*(l%LPW_F))) & (MASK_F<<(_Af*(l%LPW_F)));
    break;
  case 4:
    c[l/LPW_N] &= ~(MASK_N<<(_An*(l%LPW_N)));
    c[l/LPW_N] |= (v<<(_An*(l%LPW_N))) & (MASK_N<<(_An*(l%LPW_N)));
    break;
  }
}


// this function calculates the time spent in each patch by an individual
// time spent in a patch is calculated as a ratio of preference for the resource present in the current patch
// and to the sum pf preferences to the other patches
// with the _GENEDRIVE OPTION, it is equal time spend in all the available patches within their roaming range
static void TimeSpentRoam(pIndividual i, pPatch p)
{
  int j;
  for(j=0; j < (p->ddne); j++){
    p->dne[j]->Other[p->dne[j]->nothers].ind = i;  // this individual is saved as 'other' in another patch
    p->dne[j]->nothers++;                          // increments the no of other individuals in the patch
  }
}
