#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sched.h>
#include <sys/ipc.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include "parameters.h"
#include "random.h"
#include "niche.h"
#if _ARBITER
#include "arbiter.h"
#endif

#if _HABITAT
extern u64b     Habitat_logs;
extern FILE    *Habitatf;
#endif

#if _DISTANCE
extern u64b     Distance_logs;
extern FILE    *Distancef;
int             distancedata[12][_DRANGE+1];
#endif

int             InocGen[] = _InocGen;
int             InocCoor[_InocPatchNo][2] = _InocCoor;  // converts the data to coordinate pairs
RndState        R;

// Mutex and condition variables for thread syncronization
extern volatile pthread_cond_t  Master,Slaves;
extern volatile pthread_mutex_t Mutex;
extern volatile int             Go[_CPUS],Done[_CPUS];
extern volatile pThreadData     Threads[_CPUS];

// Random Seed, and starting space state.  Copied by the threads
extern          u64b            Seed;
extern          Space           Start;

// Array based parameters (to be coppied to threads) 
// Could these arrays be turned into constants?
//extern        double          Sigma_s[];
//extern        double          Sigma_c[];
//extern        double          a_i[];
extern          double          Xi[];
//extern        double          _LOFPROB;

// The current and total generation numbers
extern          u64b            Gen;
extern          u64b            NGen;
#if _EXPOSURE
extern          u64b            gen_exp;
#endif
extern          int             OffspringPopSize;
extern          int             Mortality;
extern          int             SinglePaternity;
extern          int             MultiplePaternity;
extern          int             MatingMotherwithDrive;
extern          int             MatingFatherwithDrive;



// Code shared with niche.c
#include "shared.h"

/***************************************************************************************
 * Usefull "utility" functions
 ***************************************************************************************/

static inline void ArbiterProgress()
{
#if _ARBITER
  arbiter_progress((((double)Gen)/NGen)*100.0);
#endif
}


/*
// in random.c; Returns choice from poisson distribution with mean mu
static int Poisson(double mu)
{

  double p = exp( -mu ), q = 1;
  int    n=-1;

  do q *= U01(&R), n++;
  while( (q > p) || (q == p) );

  return n;

}
*/

/***************************************************************************************
 * Functions dealing with the creation and syncronization of threads
 ***************************************************************************************/

// Signals that the thread td has completed and then waits for a go (or exit) signal
static int SignalDone(pThreadData td) 
{
  // Signal done and wait for Go (or exit) singal
  pthread_mutex_lock((pthread_mutex_t*)&Mutex);
  Done[td->id] = 1;
  pthread_cond_signal((pthread_cond_t*)&Master);
  while( !Go[td->id] )
    pthread_cond_wait((pthread_cond_t*)&Slaves, (pthread_mutex_t*)&Mutex);
  Go[td->id] = 0;
  pthread_mutex_unlock((pthread_mutex_t*)&Mutex);

  // For now there is no exit condition
  return 1;
}


/***************************************************************************************
 * Functions dealing directly with maintenance of Individuals
 ***************************************************************************************/

// I'm commenting out these functions as they are not used:
// this will keep the compiler warnings to a minimun.
#if 0

// Returns the number of individuals used in the "next generation"
// Note: This count _will_ include wanderers if there are any.
static u64b nNextIndividuals()
{
  u64b i,m;

  for(m=i=0; i<_CPUS; i++)
    m += Threads[i]->nINext;
  return m;
}

// Returns the number of individuals used in the entire simulation.
// Note: This count _will_ include wanderers if there are any.
static u64b nIndividuals()
{
  u64b i,m;

  for(m=i=0; i<_CPUS; i++)
    m += Threads[i]->nICurrent + Threads[i]->nINext;
  return m;
}
#endif


// Fills in traits from characteristic array
static void SyncIndividual(pThreadData td, pIndividual i)
{

  u64b j,k;

  // Sync patch related traits from characteristic array
  for(j=0;j<_k;j++) {
    for(i->x[j]=0.,k=_Le*j; k<_Le*(j+1); k++)
      i->x[j] += ((i->x0[k/LPW_E]>>(_Ae*(k%LPW_E)))&MASK_E) + ((i->x1[k/LPW_E]>>(_Ae*(k%LPW_E)))&MASK_E);
    //for(i->y[j]=0.,k=_Lp*j; k<_Lp*(j+1); k++)
      //i->y[j] += ((i->y0[k/LPW_P]>>(_Ap*(k%LPW_P)))&MASK_P) + ((i->y1[k/LPW_P]>>(_Ap*(k%LPW_P)))&MASK_P);
    // Scale traits
    i->x[j] /= MASK_E*(_Le<<1);
    //i->y[j] /= MASK_P*(_Lp<<1);
  }

  // Sync and scale maker
  for(i->k=0.,k=0; k<_Lk; k++)
    i->k += ((i->k0[k/LPW_K]>>(_Ak*(k%LPW_K)))&MASK_K) + ((i->k1[k/LPW_K]>>(_Ak*(k%LPW_K)))&MASK_K);
  i->k /= MASK_K*(_Lk<<1);
  
  // Sync and scale mating like in others [0,1]
  for(i->m=0.,k=0; k<_Lm; k++)
    i->m += ((i->m0[k/LPW_M]>>(_Am*(k%LPW_M)))&MASK_M) + ((i->m1[k/LPW_M]>>(_Am*(k%LPW_M)))&MASK_M);
  i->m /= MASK_M*(_Lm<<1);

  /*
  // old version
  // Sync and scale mating (in interval [-1,1])
  for(i->m=0.,k=0; k<_Lm; k++)
    i->m += ((i->m0[k/LPW_M]>>(_Am*(k%LPW_M)))&MASK_M) + ((i->m1[k/LPW_M]>>(_Am*(k%LPW_M)))&MASK_M);
  i->m = (i->m/(MASK_M*_Lm)) - 1.0;
  */

  // Sync and scale mating preference
  for(i->f=0.,k=0; k<_Lf; k++)
    i->f += ((i->f0[k/LPW_F]>>(_Af*(k%LPW_F)))&MASK_F) + ((i->f1[k/LPW_F]>>(_Af*(k%LPW_F)))&MASK_F);
  i->f /= MASK_F*(_Lf<<1);

    
#if _SEXCHROMOSOME
  // Sync and scale sex chromosomes (in interval [0,1], and ONLY 0 and 1)
  for(i->s=0,k=0; k<_Ln; k++)
    i->s += ((i->z0[k/LPW_N]>>(_An*(k%LPW_N)))&MASK_N) + ((i->z1[k/LPW_N]>>(_An*(k%LPW_N)))&MASK_N);
    i->s = (i->s/(MASK_N*_Ln)) - 0.5;
#endif
}

// Installs individual i into patch p in a safe manner
static int InstallIndividual(pPatch p, pIndividual i)
{
  if(p->ni == _MAX_INDIVIDUALS) return 0;
  p->i[p->ni++] = i;            return 1;
}

// Returns a pointer to an available individual structure.  NULL for none.
static pIndividual GetNextIndividual(pThreadData td)
{
  // Record where this individual was allocated from
  td->INext[td->nINext]->id = td->nINext;
  return td->INext[td->nINext++];
}

// Releases an Individual back to the its individual array
static void ReleaseNextIndividual(pThreadData td, pIndividual i)
{
  // Swap last in use with i and decriment
  td->INext[i->id]      = td->INext[--td->nINext];
  td->INext[td->nINext] = i;

  // Update allocated from (id) field
  td->INext[i->id]->id   = i->id;
}

// Clears all individuals by restoring the initial sate of individual allocation
static void ReleaseAllCurrentIndividuals(pThreadData td)
{
  td->nICurrent = 0;
  memcpy(td->ICurrent, td->ICurrentCleared, MAX_THREAD_POP*sizeof(pIndividual));

  // this is to shuffle the ind in the memory (to prevent potential first-come-first-serve bias)
  {
    int         i;
    pIndividual ind;

    for(i=0; i<MAX_THREAD_POP; i++) {
      ind = td->ICurrent[(MAX_THREAD_POP-1)-i];
      td->ICurrent[(MAX_THREAD_POP-1)-i] = td->ICurrent[i];
      td->ICurrent[i] = ind;
    }
  }
}

// I'm commenting out these functions as they are not used:
// this will keep the compiler warnings to a minimun.
#if 0
// Clears all individuals by restoring the initial sate of individual allocation
static void ReleaseAllNextIndividuals(pThreadData td)
{
  td->nINext = 0;
  memcpy(td->INext, td->INextCleared, MAX_THREAD_POP*sizeof(pIndividual));
}
#endif


// Copies all data in a safe manner from Individual s to d
static void CopyIndividual(pIndividual d, pIndividual s)
{
  u64b did;

  // Save id field
  did  = d->id;

  // Copy data
  memcpy(d,s,sizeof(Individual));

  // Re-link and restore id field
  LinkIndividual(d);
  d->id  = did;
}

/***************************************************************************************
 * Functions which compute probabilities, properties, or preferences of individuals
 ***************************************************************************************/

// effective population size, the number of individuals present in roaming range (for fertility)
static double Ne(pThreadData td, pIndividual i, pPatch p)
{
  int     j;
  double  c = 0.;

  
  // Go through the range of the individual
  /*
  for(j=0; j<(p->rrne); j++) {
    c += p->rne[j]->ni;
  }*/
  
  // Transform p (which is a pointer into next) into "current" version of p
  p = &(td->Current[p->lat][p->lon]);
  
  for(j=0; j<(p->ddne); j++) {
    c += (p->dne[j]->ni);
    printf("(%d, %d) %d, ",p->dne[j]->lat, p->dne[j]->lon,  p->dne[j]->ni);
  }
  printf("\n");
  printf("(%d, %d): %lf (range: %d)\n", p->lat, p->lon, c, p->ddne);
  return c;
}

// effective population size, the number of individuals present in dispersal range (for dispersal)
static double NeDispersal(pThreadData td, pPatch p)
{
  int     j;
  double  c = 0.;

  // Transform p (which is a pointer into next) into "current" version of p
  p = &(td->Current[p->lat][p->lon]);
    
  // Go through the range of the individual
  for(j=0; j<(p->ddne); j++) {
    // loop through and add all the number of individuals in the roaming range
    c += p->dne[j]->ni;
  }
  return c;
}


// The number of offspring to be produced, this is the seed for Poisson distribution
static double w(pThreadData td, pIndividual i, pPatch p)
{
  double r, Kz;

  // with gene-drive option, the function K(), which is actually the fitness, is turned off
  // everybody has the same fitness, i.e. 1 (see shared.h)
  
  // donut model
  //if( (p->lat <= 9 || p->lat >= 54 || p->lon <= 9 || p->lon >= 54) ){
  //  Kz = 3*_K0;
  //} else {
  //  Kz = _K0;
  //}
  
  // Transform p (which is a pointer into next) into "current" version of p
  p = &(td->Current[p->lat][p->lon]);
  
  Kz = _K0;
  
  /*
  // in the gene-drive model, Ne() counts all the individuals in roaming range (i.e. no timespent)
  // so K is adjusted with multiple patches
  // old version where r is calculated at the entire range
  double rr;
  rr = ( ((double)(Ne(td,i,p))) / ((double)(Kz * p->ddne)) );
  printf("(%d, %d) old r is: %lf\n", p->lat, p->lon, rr);
  */
  r = ( ((double)(p->ni)) / ((double)(Kz)) );
  //printf("(%d, %d) r is: %lf with N: %d\n", p->lat, p->lon, r, p->ni);

  /*
  // print for a sanity check
  double f;
  f = (1.0/(1.0+(((double)(_b))/2.0-1.0)*r));
  printf("\n'v' in beverton-holt is: %f\n",f);
  */
  
  return (1.0/(1.0+(((double)(_b))/2.0-1.0)*r));
  //double mmax = 3.58;
  //double mm = 0.2;
  //return (mmax - (mmax - mm)*r);
}



/***************************************************************************************
 * Functions which represent the conceptual steps required for the simulation (+ helpers)
 ***************************************************************************************/

//#if SKIP_MUTATION
// Fills in mutation distributions (->MutateDist and ->MutateSkip)
// below called as: MutateProb(td, 0, PLENGTH_E, FIELDS_E, _mue, _Ae, MASK_E);
// which are typically MutateProb(td, 0, PLENGTH_E=4, FIELDS_E=32, _mue=0.00001, _Ae=1, 11111111);
// which are typically MutateProb(td, 0, PLENGTH_M=1, FIELDS_M=8,  _mum=0.00001, _Am=1, 11111111);
/*
static void MutateProb(pThreadData td, int w, int plength, int fields, double mu, int abits, u64b mask)
{
  int i,j,x,e,d;  double s,sk;

  // If no mutation, just skip all this mess
  if(!mu) return;

  // Setup ->MutateDist
  // Go through all possible 8-bit source patterns
  for(x=0; x<256; x++) {
    // Allocate space for 256 destinations and initialize their probabilities to 0
    td->MutateDist[w][x] = allocdist(256);
    memset(td->MutateDist[w][x]->p, 0, 256*sizeof(double));
    // Traverse all potential mutation results
    
    for(i=0,s=0.; i<256; i++) {
      // Consider all loci within this byte
      for (j=e=0; j < fields; j++){
        // If this destination is a non-neighbor, it is not a valid mutation destination.
        // Continue to next destination byte.
        if ((d = abs(((i>>(abits*j))&mask)-((x>>(abits*j))&mask))) > 1)
            goto next_byte;
        // Valid destination, so increment the mutation event counter if d is 1
        e += d;
      }
      if (e) s += (td->MutateDist[w][x]->p[i] = pow(mu,e)*pow(1.-mu,fields-e));
      next_byte:
      continue;
    }
    // Call initdist on this distribution
    initdist(td->MutateDist[w][x],s);
  }

  // Setup ->MutateSkip
  // Allocate space for skip distribution: 
  // skip value of {0...,PLENGTH} bytes, where PLENGTH == no mutation
  td->MutateSkip[w] = allocdist(plength+1);
  // Initialize sk to the probability that a byte will be skipped
  sk = pow(1.-mu, (double)fields);
    
  // Consider all skip possibilities which will result in mutation and fill in the appropriate probability
  for(i=0,s=0.; i<plength; i++)
    s += td->MutateSkip[w]->p[i] = pow(sk,(double)i) * (1.-sk);
        
  td->MutateSkip[w]->p[plength] = 1.-s;
  // Call initdist on this distribution
  initdist(td->MutateSkip[w],1.);
}
*/
/*
#else

// Fills in mutation distribution ->MutateDist
// below called as: MutateProb(td, 0, PLENGTH_E, FIELDS_E, _mue, _Ae, MASK_E);
static void MutateProb(pThreadData td, int w, int plength, int fields, double mu, int abits, u64b mask)
{
  int i,j,x,d;  double s,f,h,nf,nh;

  // If no mutation, just skip all this mess
  if(!mu) return;

  // Setup ->MutateDist
  // Go through all possible 8-bit source patterns
  for(x=0; x<256; x++) {
    // Allocate space for 256 destinations and initialize their probabilities to 0
    td->MutateDist[w][x] = allocdist(256);
    memset(td->MutateDist[w][x]->p, 0, 256*sizeof(double));
    for(s=0.,j=0; j<256; j++) {                                            // traverse potential mutation result
      for (i = f = h = nf = nh = 0.; i < fields; i++){                     // consider the multi-bit loci wihtin a byte
	if ((d = abs(((j>>(abits*i))&mask)-((x>>(abits*i))&mask))) > 1)        // if mutate to non-neighbor...
	  goto e;                                                              // invalid mutation result
	if (abits > 1){                                                        // some alleles can mutate both directions
	  if ( (!((x>>(abits*i))&mask)) || (((x>>(abits*i))&mask)==mask) ) {   // if at boundary...
	    if (d) h++;                                                        // count mutation event
	    else   nh++;                                                       // count non-mutation event
	  } else {                                                             // not at boundary...
	    if (d) f++;                                                        // count mutation event
	    else   nf++;                                                       // count non-mutation event
	  }
	} else {                                                               // can only mutate one direction
	  f += d;                                                              // count mutation events
	}
      }
      if (abits> 1) s += (td->MutateDist[w][x]->p[j] = pow(mu,f)*pow(1.-mu,nf)*pow(mu/2.,h)*pow(1.-mu/2.,nh));
      else          s += (td->MutateDist[w][x]->p[j] = pow(mu,f)*pow(1.-mu,fields-f));
    e:                                                                     // consider next byte value
      continue;
    }
    initdist(td->MutateDist[w][x],s);
  }
}
#endif
*/

/*
  Initializes neighbor mutation
*/
/*
static void MutateInit(pThreadData td)
{
  // Setup skip and mutation distributions for ->e and ->p
  MutateProb(td, 0, PLENGTH_E, FIELDS_E, _mue, _Ae, MASK_E);
  MutateProb(td, 1, PLENGTH_P, FIELDS_P, _mup, _Ap, MASK_P);

  // Setup skip and mutation distributions for ->m and ->k
  MutateProb(td, 2, PLENGTH_M, FIELDS_M, _mum, _Am, MASK_M);
  MutateProb(td, 3, PLENGTH_K, FIELDS_K, _muk, _Ak, MASK_K);
  // Setup skip and mutation distributions for ->f
  MutateProb(td, 4, PLENGTH_F, FIELDS_F, _muf, _Af, MASK_F);

  // Setup skip and mutation distributions for ->z
  MutateProb(td, 5, PLENGTH_N, FIELDS_N, _mun, _An, MASK_N);

}
*/
/*
  Returns the number of flips (0 <--> 1) in crossover mask x; clone event <--> x = 255
  Used in CrossProb with argument i=0; i<256
  byte=unsigned char;
*/
static double CrossFlips(byte x) 
{
  /* initial 0 <--> previous cross point */
  int a = 1-(x&1), i = 7; 

  for ( ; i--; x >>= 1) a += (x&1)^((x>>1)&1);
  return (double) a;
}

/*
  Returns 1 if crossover mask x is valid (does not split a multi-bit loci), 0 on invalid
*/
static int CrossValid(byte x, int abits)
{
  int a,i;

  /* Move through all bits */
  for (i=a=0; x; x>>=1){
    /* Sum the number of ones */
    a += x&1;
    /* Check to see if we are done with a locus here */
    if (!((++i)%abits)){
      /* If this locus is split, return 0 */
      if (a%abits) 
	return 0;
      /* Reset a and keep searching */
      a = 0;
    }
  }

  /* Test last locus and return valid/not valid */
  return ((a%abits)? 0: 1);
}

// Fills in crosover distribution
static void CrossProb(pThreadData td, int w, int fields, int abits)
{
  int i;  double s;

  /* Allocate space for all possible crossover masks */
  td->CrossDist[w] = allocdist(256);
  // note: The C library function void *memset(void *str, int c, size_t n)
  // copies the character c (an unsigned char) to the first n characters of the string pointed to,
  // by the argument str.
  // sizeof(double) = 8 bytes
  // so this is making everything zero?
  memset(td->CrossDist[w]->p, 0, 256*sizeof(double));

  /* Consider all crossover masks */
  for (s=0.,i=0; i<256; i++){
    /* Make sure mask is valid and fill in it's prob */
    // when is it not valid???
    if (CrossValid(i,abits))
      s += td->CrossDist[w]->p[i] = pow(Xi[w],CrossFlips(i)) * pow(1.-Xi[w],fields-CrossFlips(i));
    //printf("w, Xi[w], i, s, td->CrossDist[w]->p[i]: %i, %lf, %i, %0.18llf, %0.18llf\n",w, Xi[w], i,s, td->CrossDist[w]->p[i]);
    }
    /* Call initdist to initialize this distribution */
    /*
     in random.c it says: "Note: d->p must have d->n elements which sum to s on entry to initdist.
    The elements of d->p and d->a are overwritten by the initialization process."
    */
  initdist(td->CrossDist[w], s);
}

// Initialize crossover
static void CrossInit(pThreadData td)
{
  // Initialize distributions ->e and ->p
  
  CrossProb(td, 0, FIELDS_E, _Ae);
  // Initialize distributions ->m and ->k
  CrossProb(td, 1, FIELDS_M, _Am);
  CrossProb(td, 2, FIELDS_K, _Ak);
  // Initialize distributions ->f
  CrossProb(td, 3, FIELDS_F, _Af);
  // Initialize distributions ->z
  CrossProb(td, 4, FIELDS_N, _An);
}


// Returns a pointer to the patch ox,oy spaces from x,y, or NULL if out of bounds
// This is not a generic function. It is mostly specific to NeighboringDispersePatches()
static pPatch ValidNeighborPatch(pThreadData td, u64b x, u64b y, int ox, int oy)
{

  // Width bounds (Only bound first left and last right)
  if( (x == 0)        && (ox < 0) )     return NULL;
  if( (x == WIDTH-1)  && (ox > 0) )     return NULL;
  if( ((x + ox) < 0) )                  return NULL;
  if( ((x + ox) > (WIDTH-1)) )          return NULL;

  // Height bounds
  if( (y == 0)        && (oy < 0) )     return NULL;
  if( (y == HEIGHT-1) && (oy > 0) )     return NULL;
  if( ((y + oy) < 0) )                  return NULL;
  if( ((y + oy) > (HEIGHT - 1)) )       return NULL;
  
  return &td->Next[x+ox][y+oy];
}


// Returns a pointer to the patch ox,oy spaces from x,y, or NULL if out of bounds
// This is not a generic function. It is mostly specific to NeighboringDispersePatches()
static int ValidMatePatch(pThreadData td, u64b x, u64b y, int ox, int oy)
{

  // Width bounds (Only bound first left and last right)
  if( (x == 0)        && (ox < 0) )     return 0;
  if( (x == WIDTH-1)  && (ox > 0) )     return 0;
  if( ((x + ox) < 0) )                  return 0;
  if( ((x + ox) > (WIDTH-1)) )          return 0;

  // Height bounds
  if( (y == 0)        && (oy < 0) )     return 0;
  if( (y == HEIGHT-1) && (oy > 0) )     return 0;
  if( ((y + oy) < 0) )                  return 0;
  if( ((y + oy) > (HEIGHT - 1)) )       return 0;

  return 1;
}


// Fills in the neighboring dispersal patches field of td->Next[x][y].
static void NeighboringDispersePatches(pThreadData td, u64b x, u64b y)
{

  int i,j;
  int dist=0;
  int count=0;
  
  for(td->Next[x][y].ddne=0,i=-((_DRANGE - 1.) / 2.); i<((_DRANGE - 1.) / 2.)+1.; i++){
    for(j=-((_DRANGE - 1.) / 2.); j<((_DRANGE - 1.) / 2.)+1.; j++) {
      if( (td->Next[x][y].dne[td->Next[x][y].ddne] = ValidNeighborPatch(td,x,y,i,j)) ){
        // density and distance dependent based dispersal
        // rather complicated way to calculate probability but basically normalized
        // by the number of patches, say if 3x3, dis=1, there are 8 patches with dis=1
        // but if you repeat patches with dis=1 eight times, the probability of
        // picking it is higher compared to the centre, which is only 1 patch.
        // so the total of all eight patches should have a prob=1+1. where center patch prob = 1
#if _DISTANCEDENSITYDISPERSAL
        // find distance
        if( abs(i) > abs(j) ){
          dist = abs(i);
        } else {
          dist = abs(j);
        }
      
        //double denom = 0.;
        double range = 0.;
        range = (((double)(_DRANGE)) - 1.)/2. + 1;
        
        // to prevent divison by zero error
        // *****************************************but it is  not correct!!!
        //if( dist == 0 ){
        //  denom = 1.;
        //} else {
        //  denom = ((double)(8*dist));
        //}
      
        td->Next[x][y].disdist[count] = dist; // distance data
        //td->Next[x][y].lowdp[count]   = ( ((double)(dist+1))*((double)(dist+1))*((double)(dist+1))/denom ); // low D probability
        //td->Next[x][y].highdp[count]  = (range - dist)*(range - dist)*(range - dist)/denom; // high D probability
        //printf("(%d, %d) %d, %lf, %lf\n",x,y,td->Next[x][y].disdist[count],td->Next[x][y].lowdp[count],td->Next[x][y].highdp[count]);
        count++;
#endif
        td->Next[x][y].ddne++;
      }
    }
  }
  //printf("\nno of patches to disperse to: %d", td->Next[x][y].ddne);
}


// Swaps Current and Next generations and deps
static void SwapGenerations(pThreadData td)
{

  u64b t;  pIndividual *tp;  Space ts;  

  // Swap spaces
  ts = td->Current;
  td->Current = td->Next;
  td->Next = ts;

  // Swap Individual resources
  t = td->nICurrent;
  td->nICurrent = td->nINext;
  td->nINext = t;

  tp = td->ICurrentCleared;
  td->ICurrentCleared = td->INextCleared;
  td->INextCleared = tp;

  tp = td->ICurrent;
  td->ICurrent = td->INext;
  td->INext = tp;
}


// Initializes any data structures needed by a thread (mostly [if not all] local)
static void InitThread(pThreadData td)
//void InitThread(pThreadData td)
{
  u64b i,j,k,l;
  
  // Make a copy of array based parameters
  //ALLOC(td->Sigma_s,     _k*sizeof(double), "Could not allocate Sigma_s!\n");
  //ALLOC(td->Sigma_c,     _k*sizeof(double), "Could not allocate Sigma_c!\n");
  //ALLOC(td->a_i,         _k*sizeof(double), "Could not allocate a_i!\n");
  //ALLOC(td->log_Sigma_s, _k*sizeof(double), "Could not allocate log_Sigma_s!\n");
  //memcpy(td->Sigma_s,Sigma_s,_k*sizeof(double));
  //memcpy(td->Sigma_c,Sigma_c,_k*sizeof(double));
  //memcpy(td->a_i,    a_i,    _k*sizeof(double));

  // Precompute the log of Sigma_s
  //for(i=0; i<_k; i++)
  //  td->log_Sigma_s[i] = log(td->Sigma_s[i]);

  // Initialize Individuals and set up the "cleared" states
  ZALLOC(td->ICurrentData,   MAX_THREAD_POP*sizeof(Individual),  "Could not allocate IndividualData! (Cur)\n");
  ALLOC(td->ICurrent,        MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Individuals! (Cur)\n");
  ZALLOC(td->INextData,      MAX_THREAD_POP*sizeof(Individual),  "Could not allocate IndividualData! (Next)\n");
  ALLOC(td->INext,           MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Individuals! (Next)\n");
  ALLOC(td->ICurrentCleared, MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Cleared State! (Cur)\n");
  ALLOC(td->INextCleared,    MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Cleared State! (Next)\n");
  for(td->nICurrent=td->nINext=i=0; i<MAX_THREAD_POP; i++){
    td->ICurrent[i]      = td->ICurrentData+i;
    td->INext[i]         = td->INextData+i;
    LinkIndividual(td->ICurrent[i]);
    LinkIndividual(td->INext[i]);
  }
  memcpy(td->ICurrentCleared, td->ICurrent, MAX_THREAD_POP*sizeof(pIndividual));
  memcpy(td->INextCleared,    td->INext,    MAX_THREAD_POP*sizeof(pIndividual));

  // Malloc for our patch spaces
  ALLOC(td->Current,(WIDTH+2)*sizeof(Patch*),   "Could not allocate Current patch space!\n");
  for(i=0; i<(WIDTH+2); i++)
    ZALLOC(td->Current[i],HEIGHT*sizeof(Patch), "Could not allocate Current patch space!\n");
  ALLOC(td->Next,(WIDTH+2)*sizeof(Patch*),      "Could not allocate Next patch space!\n");
  for(i=0; i<(WIDTH+2); i++)
    ZALLOC(td->Next[i],HEIGHT*sizeof(Patch),    "Could not allocate Next patch space!\n");

  // Initialize the random number generator
  initrand(&R, Seed);

  
  // Figure out what patches to copy from start
#if _CPUS == 1
  // Only one thread:    Do not copy any wandering zones
  k = 1, l = WIDTH+1;
#else
  switch(td->id){
  case 0:         // Left-most thread:   Only copy right wandering zone
    k = 1, l = WIDTH+2;
    break;
  case (_CPUS-1): // Right-most thread:  Only  copy left wandering zone
    k = 0, l = WIDTH+1;
    break;
  default:        // Middle threads:     Copy both wandering zones
    k = 0, l = WIDTH+2;
  }
#endif

  // Actually copy them
  for(i=0; i<WIDTH; i++){
    // Copy into current
    memcpy(td->Current[i],    Start[i],          HEIGHT*sizeof(Patch));
    // Next matches current (no individuals)
    memcpy(td->Next[i],       td->Current[i],    HEIGHT*sizeof(Patch));
    for(j=0; j<HEIGHT; j++) td->Next[i][j].ni = 0;
  }

  // Copy extinction information
  /*
  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++) {
      td->Current[i][j].en = Start[i][j].en;
      ALLOC(td->Current[i][j].ei, td->Current[i][j].en*sizeof(ExtinctionInfo),
            "Could not allocate space for extinction information! (Wt)\n");
      memcpy(td->Current[i][j].ei, Start[i][j].ei,
             td->Current[i][j].en*sizeof(ExtinctionInfo));
      // Just make next point to current's extinction events
      td->Next[i][j].ei = td->Current[i][j].ei;
      td->Next[i][j].en = td->Current[i][j].en;
    }
  }*/

  // Fill in patch individuals
  for(i=0; i<WIDTH; i++){
    // *for(i=1; i<(WIDTH+1); i++) {
    // Copy individuals from start
    for(j=0; j<HEIGHT; j++){
      for(k=0; k<td->Current[i][j].ni; k++){
        //printf("\nenters the loop. start: %d, current: %d", Start[i][j].ni, td->Current[i][j].ni);
        CopyIndividual(td->Current[i][j].i[k]=GetCurrentIndividual(td),Start[i][j].i[k]);
        // Start[][]'s individuals are unsynced
        SyncIndividual(td,td->Current[i][j].i[k]);
      }
    }
  }

  // Now that patches are copied, pre-compute neighbors, but Since NeighboringPatches()
  // only fills in the Next generation, swap and process twice.
  SwapGenerations(td);
  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++){
      NeighboringDispersePatches(td,i,j);
    }
  }

  SwapGenerations(td);
  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++){
      NeighboringDispersePatches(td,i,j);
    }
  }

  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++) {
      td->Current[i][j].nothers=0;
    }
  }

  // calculates who is using which patch for how long
  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++) {
      for(k=0,l=td->Current[i][j].ni; k < l; k++) {
        TimeSpentRoam(td->Current[i][j].i[k], &td->Current[i][j]);
        //TimeSpentMate(td->Current[i][j].i[k], &td->Current[i][j]);
      }
    }
  }

  // Setup recombination and mutation distributions
  //MutateInit(td);
  CrossInit(td);
}

#if ENDIAN == NICHE_BIG_ENDIAN
// This macro will do an index translation that should be 
// functionally equvalent to byte-swapping.  
static int BSWAP[32] = {
 7,  6,  5,  4,  3,  2,  1,  0,
 15, 14, 13, 12, 11, 10, 9,  8,
 23, 22, 21, 20, 19, 18, 17, 16,
 31, 30, 29, 28, 27, 26, 25, 24
};
#define Endian(x) (BSWAP[x])
#else
#if ENDIAN == NICHE_LITTLE_ENDIAN
#define Endian(x) (x)
#else
#error ENDIAN _must_ be set to either NICHE_BIG_ENDIAN or NICHE_LITTLE_ENDIAN!
#endif
#endif

// Handles recombination and mutation for a single destinaation chromosome
static void RecombineEco(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)){ t=x; x=y; y=t;}

#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_E; i++){
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[0]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_E; i++){
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[0]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
//gene-drive: #endif
  }
#endif
}


// Handles recombination and mutation for a single destinaation chromosome
static void RecombineNeutral(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }
  
  //printf("before - z, x, y: %d, %d, %d\n", *z, *x, *y);
  
#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_N; i++){
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[4]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_N; i++){
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[4]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
//gene-drive: #endif
    //printf("after - z, x, y: %d, %d, %d\n", *z, *x, *y);
  }
#endif
}


// Handles recombination and mutation for a single destinaation chromosome
static void RecombineMating(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }
  //printf("BEFORE x, y, z: %llu, %llu, %llu\n", *x, *y, *z);
    
#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_M; i++){
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[1]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_M; i++){
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[1]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
    // if( x[Endian(i)] == 255) printf("x: %llu\n", x[Endian(i)]);
    // printf("x, y, z: %llu, %llu, %llu\n", x[Endian(i)], y[Endian(i)], z[Endian(i)]);
//gene-drive: #endif
  }
#endif
}


// Handles recombination and mutation for a single destinaation chromosome
static void RecombineMarker(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }
    
#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_K; i++) {
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[2]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_K; i++) {
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[2]);
    //printf("r: %llu/n", r);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
// gene-drive: #endif
  }
#endif
}


// Handles recombination and mutation for a single destination chromosome
static void RecombineFPref(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }

#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_F; i++) {
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[3]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }

#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_F; i++) {
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[3]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
//gene-drive: #endif
  }
#endif
}


static double GetFemaleGermLine()
{
  distribution *gdist;
  double germ=4.;
  double probsum = 0.0;
  gdist = allocdist(4);
  gdist->n = 4;
  

  // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
  probsum += ( gdist->p[0] = ((double)(0.5 + _FEMALECUTTINGPROB/2.0 - (_FEMALECUTTINGPROB * _FEMALENHEJPROB)/2.0)) );// drive (allele 2, or 5 after exp/homer-hs)
  probsum += ( gdist->p[1] = ((double)(0.5 * _FEMALECUTTINGPROB * _FEMALENHEJPROB * _FEMALELOFPROB)) );            // nhej-LOF (allele 4)
  probsum += ( gdist->p[2] = ((double)((1. - _FEMALECUTTINGPROB)/2.0)) );                                          // drive fails: wildtype (allele 1)
  probsum += ( gdist->p[3] = ((double)(0.5 * _FEMALECUTTINGPROB * _FEMALENHEJPROB * (1. - _FEMALELOFPROB))) );     // nhej-functional (allele 3)
  
  initdist(gdist,probsum);
  germ = drand(&R, gdist);
  
  //if( germ == 3.0){ printf("inside the function germ: %lf\n", germ); }
  freedist(gdist);
  return germ;
}

static double GetMaleGermLine()
{
  distribution *gdist;
  double germ=4.;
  double probsum = 0.0;
  gdist = allocdist(4);
  
  // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
  probsum += ( gdist->p[0] = (0.5 + _MALECUTTINGPROB/2 - (_MALECUTTINGPROB * _MALENHEJPROB)/2.0) );// drive (allele 2)
  probsum += ( gdist->p[1] = (0.5 * _MALECUTTINGPROB * _MALENHEJPROB * _MALELOFPROB) );            // nhej-LOF (allele 4)
  probsum += ( gdist->p[2] = ((1. - _MALECUTTINGPROB)/2.0) );                                      // drive fails: wildtype (allele 1)
  probsum += ( gdist->p[3] = (0.5 * _MALECUTTINGPROB * _MALENHEJPROB * (1. - _MALELOFPROB)) );     // nhej-functional (allele 3)

  
  gdist->n = 4;
  initdist(gdist,probsum);
  germ = drand(&R, gdist);
  //if( germ == 3.0){ printf("inside the function germ: %lf\n", germ); }
    //if( germ == 2.0){ printf("m germ:%lf- %lf, %lf, %lf, %lf\n", germ, gdist->p[0], gdist->p[1], gdist->p[2], gdist->p[3]); }
  //printf("inside the function germ: %lf\n", germ);
  freedist(gdist);
  return germ;
}


// Combines parents m and f and fills in child c with correct recombination and mutation
static void FastRecombine(pThreadData td, pIndividual c, pIndividual m, pIndividual f)
{
  
  ArbiterProgress();
  
  // Random parent swap (i am silencing this because order is useful now)
  // pIndividual t;
  //if (rnd(&R, 2)) { t = m; m = f; f = t; }
  //if (m->s != 0) printf("there is an error in FastRecombine(), m should be a male\n");
  
  // Even though these are called Recombine, it only deals with
  // chromosomal recombination shuffling loci on the same 'chromosome' i.e. trait,
  // independent traits that are on different chromosomes, i.e. are independent anyway
  // there is no cross-over bwtween loci
  // as _Xi is defined as 0 in parameters.h
  // could be turned on, but NOT in drive related loci
  RecombineEco     (td, (byte*)(c->x0), (byte*)(m->x0), (byte*)(m->x1));
  RecombineEco     (td, (byte*)(c->x1), (byte*)(f->x0), (byte*)(f->x1));
  RecombineMarker  (td, (byte*)(c->k0), (byte*)(m->k0), (byte*)(m->k1));
  RecombineMarker  (td, (byte*)(c->k1), (byte*)(f->k0), (byte*)(f->k1));
  RecombineMating  (td, (byte*)(c->m0), (byte*)(m->m0), (byte*)(m->m1));
  RecombineMating  (td, (byte*)(c->m1), (byte*)(f->m0), (byte*)(f->m1));
  RecombineNeutral (td, (byte*)(c->z0), (byte*)(m->z0), (byte*)(m->z1));
  RecombineNeutral (td, (byte*)(c->z1), (byte*)(f->z0), (byte*)(f->z1));
  RecombineFPref   (td, (byte*)(c->f0), (byte*)(m->f0), (byte*)(m->f1));
  RecombineFPref   (td, (byte*)(c->f1), (byte*)(f->f0), (byte*)(f->f1));
}


#if _DISTANCE
void RecordDistance(double density, int distance)
{
  
  if( (density >= 0.0) && (density < 0.1) ) { distancedata[0][distance]++; }
  if( (density >= 0.1) && (density < 0.2) ) { distancedata[1][distance]++; }
  if( (density >= 0.2) && (density < 0.3) ) { distancedata[2][distance]++; }
  if( (density >= 0.3) && (density < 0.4) ) { distancedata[3][distance]++; }
  if( (density >= 0.4) && (density < 0.5) ) { distancedata[4][distance]++; }
  if( (density >= 0.5) && (density < 0.6) ) { distancedata[5][distance]++; }
  if( (density >= 0.6) && (density < 0.7) ) { distancedata[6][distance]++; }
  if( (density >= 0.7) && (density < 0.8) ) { distancedata[7][distance]++; }
  if( (density >= 0.8) && (density < 0.9) ) { distancedata[8][distance]++; }
  if( (density >= 0.9) && (density < 1.0) ) { distancedata[9][distance]++; }
  if( (density >= 1.0) && (density < 1.1) ) { distancedata[10][distance]++; }
  if( (density >= 1.1) )                    { distancedata[11][distance]++; }

}


void LogDistance()
{

  int     x,y;
  double  maxdis = 0.;
  
  maxdis = ( ((double)(_DRANGE)) - 1. )/2. + 1.;
  
  // Gen+1 is for viewing
  // otherwise since dispersal occurs at gen 0, it is viewed at gen 0
  // we want to see it 'after' it happened at gen 1
  //printf("Logging distance at gen: %d as %d\n", Gen, Gen+1);
  if( !(Gen%_DISTANCE) )
    fprintf(Distancef, "gen: %llu ", Gen+1);
  for( x=0; x<12; x++){
    //i-- > 0 ;
    //for( y=0; y<((int)(maxdis)); y++){
    for( y= ((int)(maxdis)); y-- >0;){
      fprintf(Distancef, "%d ", distancedata[x][y]);
      //printf("%d ", distancedata[x][y]);
    }
  }
  fprintf(Distancef, "\n");
  fflush(Distancef);
 // Increment counter
 Distance_logs++;
}
#endif


// dispersal is done w.r.t. one of the neighboring patches of the mother's
// this could be modified to make the dispersal from the patch where mating took place
static void Disperse(pThreadData td, pIndividual i, pPatch p)
{
  int     distance, k;
  static  distribution *d  = NULL;
  static  distribution *cd = NULL;
  double  s, sd, maxdis, x, normalizeddis;
  
  x             = 0.;
  normalizeddis = 0.;
  maxdis  = ((double)((_DRANGE - 1.) / 2.));
  
  ArbiterProgress();
  
  /*
  {
    u32b    rtab[55];
    int     rndx,nflg;
    double  nv;
  } RndState;
*/
  // Only allocate the neighbor distribution the first time
  /*
   //old version
  if( !d ){
    d = allocdist(_DRANGE * _DRANGE);
  }
  d->n = p->ddne;
  */
  
#if _DISTANCEDENSITYDISPERSAL
  
  d = allocdist(((int)(maxdis))+1);
  d->n = ((int)(maxdis))+1;
 
  cd = allocdist(p->ddne);
  cd->n = p->ddne;
  
  // Transform p (which is a pointer into next) into "current" version of p to get density
  p = &(td->Current[p->lat][p->lon]);
  //printf("current patch is: %d, %d with %d ind\n", p->lat, p->lon, p->ni);
  
  // density and distance dependent based dispersal
  for( s=0., k=0; k<((int)(maxdis))+1; k++ ){
    d->p[k] = 0.0;
    //printf("%lf\n", d->p[k]);
    normalizeddis = ((double)(k))/maxdis;
    x = (_DISPERSALCOEF*normalizeddis - _DENSITYCOEF*p->density);
    s += ( d->p[k] = Exp(x*x) );
    //printf("k: %d, maxdis: %d, ndis: %lf, x: %lf, Exp(x^2): %lf, s:%lf\n",k,((int)(maxdis)), normalizeddis, x, Exp(x*x), s);
    //printf("---distance: %d x: %lf, Exp(x*x): %lf, s: %lf, density: %lf\n", k, x, Exp(x*x), s, p->density);
    //printf("%lf\n", d->p[k]);
  }

  if( s != 0. ){
    initdist(d,s);
    distance = drand(&R, d);
    //printf("picked d is: %d\n", distance);
  }
  
  for( sd=0.,k=0; k<p->ddne; k++ ){
    //printf("k: %d, sd: %lf\n", k, sd);
    if( distance == p->disdist[k] ){
      sd += ( cd->p[k] = 1. );
      //printf("distance: %d, sd: %lf, lat: %d, lon: %d \n", distance, sd, p->dne[k]->lat, p->dne[k]->lon);
    } else {
      sd += ( cd->p[k] = 0.);
      //printf("sd: %lf\n",sd);
    }
  }
  //printf("final sd: %lf\n",sd);
  
  k = 0;
  if( sd != 0. ){
    //cd->n = p->ddne;
    //printf("sd: %lf\n", sd);
    initdist(cd,sd);
    k = drand(&R, cd);
    //printf("k is: %d destination is (%d, %d)\n", k, p->dne[k]->lat, p->dne[k]->lon);
  }
  
  //printf("\nfrom (%d,%d), to (%d,%d)", p->lat, p->lon, p->dne[k]->lat, p->dne[k]->lon);
#if _DISTANCE
  //printf("recording: den: %lf, dis: %d\n", p->density, distance);
  RecordDistance(p->density, distance);
#endif
#endif

#if _RAND_DISPERSAL
  // dispersal is random to one
  // Child will randomly to one of the patches within the range of the parent (female)
  d->n = p->ddne;
  k = rnd(&R, p->ddne);
  //printf("dispersal destination %d, %d\n",p->dne[k]->lat, p->dne[k]->lon);
#endif

  /*
#if _PREF_DISPERSAL
    // preference based dispersal
    for( s=0.,k=0; k<p->ddne;k++ ){
      s += (d->p[k]=P(i,p->dne[k]));    // dispersal not affected by distance
    }
    if( !s ){
      k = rnd(&R, p->ddne);
      // printf("\ndestination is %d",k);
    } else {
      initdist(d,s);
      k = drand(&R, d);
      //printf("\ndestination is %d/%d at (%d,%d)",k, p->ddne, p->dne[k]->lat, p->dne[k]->lon);
    }
#endif
    */
    // Transform p back to "next"
    
  freedist(d);
  freedist(cd);
  p = &(td->Next[p->lat][p->lon]);
  
  // sanity print
  //printf("from (%d,%d), to (%d,%d)\n", p->lat, p->lon, p->dne[k]->lat, p->dne[k]->lon);
  //printf("now in Disperse()\n");
  if(!InstallIndividual(p->dne[k],i)) {
    ReleaseNextIndividual(td,i);
    printf("release next\n");
    Error("!! HITTING UPPER BOUND ON MAX_INDS\n");
  }
}


// mates m and f, and there is fertility selection
static void Mate(pThreadData td, pIndividual f, pIndividual m, pPatch p)
{
  ArbiterProgress();
  pIndividual i;
  int         j, rn;
  double      mknockout, fknockout, nc = 0.;
  u64b        l;
  
  if ( m->s == 1 ) printf("*This ERROR shouldn't print, check Mate() m is NOT a male");
  if ( f->s == 0 ) printf("*This ERROR shouldn't print, check Mate() f is NOT a female");

  //printf("Mating btw %d(%d) and %d(%d)\n", f->id, f->s, m->id, m->s);
  /*
  // sanity print
   
  if( Gen == 210 ){
    printf("\nMating btw %d f: (%d, %d), and %d m: (%d, %d)\n", f->id, *(f->f0), *(f->f1), m->id, *(m->f0), *(m->f0));
  }
  */
    
  nc = _b*w(td,f,p);
  //nc = w(td,f,p);
  rn = Poisson(&R, nc);
  //printf("no. of offspring seed (w) for Poisson: %f - no. of offspring: %d (%d, %d)\n", nc, rn, p->lat, p->lon);
  OffspringPopSize += rn;
  
  for(j=0; j<((int)(rn)); j++) {
    // Get an individual to use as a child
    i = GetNextIndividual(td);
    i->age = 0;

    // Recombine and mutate
    FastRecombine(td, i, m, f);

    if( Gen < 174 && ((*(f->f0) == 2 && *(f->f1) == 1) || (*(f->f0) == 1 && *(f->f1) == 2)) ) MatingMotherwithDrive++;
    if( Gen < 174 && ((*(m->f0) == 2 && *(m->f1) == 1) || (*(m->f0) == 1 && *(m->f1) == 2)) ) MatingFatherwithDrive++;
      
    double fgerm = 0.;
    double mgerm = 0.;
    
    // ***************** drive
    // the following happens in the approaches
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // in _HOMING:  the drive is inserted in fertility/viability locus and is happlosufficient
    // in _HOMER:   the drive is inserted in VIABILITY locus and is happloINsufficient
    // in _HOMERKO: the drive is inserted in VIABILITY locus and is happloINsufficient
    // in _HOMERKO there is additional knockout in females (or both sexes) of a distal gene
    // if it has the drive (the allele on the other chr should be wildtype (no homing if functional resistant because already cut))
    // has to be heterozygous otherwise regular transmission
    // the following "homing" happens in all three approaches
    if( (*(f->f0) == 2 && *(f->f1) == 1) || (*(f->f0) == 1 && *(f->f1) == 2) ){
      fgerm = GetFemaleGermLine();
      //if( Gen == 122 ) { printf("female germ: %llf (%llu, %llu)\n", fgerm, *(f->f0), *(f->f1)); }
      // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
      // FEMALES f1
      if( ((int)(fgerm)) == 0 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 2, 3);
      }
      if( ((int)(fgerm)) == 1 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 4, 3);
      }
      if( ((int)(fgerm)) == 2 ){
        // never this in females if pC = 1
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 1, 3);
      }
      if( ((int)(fgerm)) == 3 ){
        // never this in females if pL = 1
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 3, 3);
      }
    }
 
    // the following happens in the approaches
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // in _HOMING:  the drive is inserted in fertility/viability locus and is happlosufficient
    // in _HOMER:   the drive is inserted in VIABILITY locus and is happloINsufficient
    // in _HOMERKO: the drive is inserted in VIABILITY locus and is happloINsufficient
    // in _HOMERKO there is additional knockout in females (or both sexes) of a distal gene
    // if it has the drive (the allele on the other chr should be wildtype (no homing if functional resistant because already cut))
    // has to be heterozygous otherwise regular transmission
    // the following "homing" happens in all three approaches
    if( (*(m->f0) == 2 && *(m->f1) == 1) || (*(m->f0) == 1 && *(m->f1) == 2) ){
      mgerm = GetMaleGermLine();
      //if( Gen == 122 ) { printf("male germ: %llf (%llu, %llu)\n", mgerm, *(m->f0), *(m->f1)); }
      // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
      // MALES f0
      if( ((int)(mgerm)) == 0 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 2, 3);
      }
      if( ((int)(mgerm)) == 1 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 4, 3);
      }
      if( ((int)(mgerm)) == 2 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 1, 3);
      }
      if( ((int)(mgerm)) == 3 ){
        // never this in males if pL = 1
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 3, 3);
      }
    }

      
#if _HOMER_HS
#if _EXPOSURE
    // the drive allele has changed from 2 to 5, which has lost the rescue, but still "homing"
    if( gen_exp > 0 && Gen >= gen_exp ){
      // get female germline
      if( (*(f->f0) == 5 && *(f->f1) == 1) || (*(f->f0) == 1 && *(f->f1) == 5) ){
        fgerm = GetFemaleGermLine();
        //if( Gen == 122 ) { printf("female germ: %llf (%llu, %llu)\n", fgerm, *(f->f0), *(f->f1)); }
        // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant, 5: drive withour rescue after exposure
        // FEMALES f1
        if( ((int)(fgerm)) == 0 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 5, 3);
        }
        if( ((int)(fgerm)) == 1 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 4, 3);
        }
        if( ((int)(fgerm)) == 2 ){
          // never this in females if pC = 1
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 1, 3);
        }
        if( ((int)(fgerm)) == 3 ){
          // never this in females if pL = 1
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 3, 3);
        }
      }
     
      // the following happens in the approaches
      // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
      // in _HOMING:  the drive is inserted in fertility/viability locus and is happlosufficient
      // in _HOMER:   the drive is inserted in VIABILITY locus and is happloINsufficient
      // in _HOMERKO: the drive is inserted in VIABILITY locus and is happloINsufficient
      // in _HOMERKO there is additional knockout in females (or both sexes) of a distal gene
      // if it has the drive (the allele on the other chr should be wildtype (no homing if functional resistant because already cut))
      // has to be heterozygous otherwise regular transmission
      // the following "homing" happens in all three approaches
      if( (*(m->f0) == 5 && *(m->f1) == 1) || (*(m->f0) == 1 && *(m->f1) == 5) ){
        mgerm = GetMaleGermLine();
        //if( Gen == 122 ) { printf("male germ: %llf (%llu, %llu)\n", mgerm, *(m->f0), *(m->f1)); }
        // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant, 5: drive withour rescue after exposure
        // MALES f0
        if( ((int)(mgerm)) == 0 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 5, 3);
        }
        if( ((int)(mgerm)) == 1 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 4, 3);
        }
        if( ((int)(mgerm)) == 2 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 1, 3);
        }
        if( ((int)(mgerm)) == 3 ){
          // never this in males if pL = 1
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 3, 3);
        }
      }
    }
#endif
#endif
      
    // ***************** knockout in distal loci
    // m locus alleles are 0: nonfunctional, 1: wildtype (_HOMERKO)
    // after the offspring is produced, knockout the genes
    // knockout happens in females in A and C; both sexes in B and D (so also include B and D in the first '#if' conditon below)
    // target for knockout is female fertility in A and B; and male fertility in C and D
    // if both A and C (or both B and D) are both on 'viability' in both sexes (Main())
#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C
    // note that B is handled by separate male knockout probability!
    // the fertility knockout is in the germline only (in ALL approach2's!!!)
    // if the mother has the drive, then the offspring will inherit '0' always!
    // ***insert probability of knockout
    // no resistance
    // check if the mother (1st cond) has the drive (2nd cond),
    // THEN knockout fertility (could also check IF it has prolactin, or just knockout anyways)
    // check if the female has the drive in 'f' loci
    if ( (f->s == 1) && ( *(f->f0) == 2 || *(f->f1) == 2 ) ){
      //knockout the 'm' locus in females (assume it inherits 'm1')
      fknockout = 0.;
      fknockout = U01(&R);
      if( fknockout <= _FEMALEKNOCKOUT ){
        for(l=0; l<_Lm; l++) SetLocus(i->m1, l, 0, 1);
      }
    }

    // the fertility knockout is in the male germline as well in B and D in addition to female knockout in A and C
    // if the father has the drive (in this additional option), then the offspring will inherit '0' always!
    // ***insert probability of knockout
    // no resistance
    // check if the father (1st cond) has the drive (2nd cond),
    // THEN knockout fertility (could also check IF it has prolactin, or just knockout anyways)
    // check if the male has the drive in 'f' loci
    if ( (m->s == 0) && ( *(m->f0) == 2 || *(m->f1) == 2 ) ){
      mknockout = 0.;
      mknockout = U01(&R);
      if( mknockout <= _MALEKNOCKOUT ){
        //knockout the 'm' locus in males (assume it inherits 'm0')
        for(l=0; l<_Lm; l++) SetLocus(i->m0, l, 0, 1);
      }
    }
#endif
    
    // ***************** viability check both HomeR and HomeRKO
#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C || _HOMER
    // // haploinsufficient: check for viability of offspring in 'f' locus
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // one copy of '4' will make it inviable (allele '2' with drive also has 'rescue')
    if ( *(i->f0) == 4 || *(i->f1) == 4  ){
      // this should work
      continue;
    }
#endif

#if _HOMERKO_C
    // // haplosufficient: check for viability of offspring in 'm' locus in C
    // alleles are 1: wildtype, 0: nonfunctional
    if ( *(i->m0) == 0 && *(i->m1) == 0  ){
      // this should work
      continue;
    }
#endif
            
#if _HOMING_C
    // // haploinsufficient: check for viability of offspring in 'f' locus
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // one copy of '4' will make it inviable (allele '2' with drive also has 'rescue')
    if ( (*(i->f0) == 2 && *(i->f1) == 2) ||
         (*(i->f0) == 2 && *(i->f1) == 4) ||
         (*(i->f0) == 4 && *(i->f1) == 2) ||
         (*(i->f0) == 4 && *(i->f1) == 4) ){
      // this should work
      continue;
    }
#endif
      
      
    // Sync, Handle dispersal and install child
    SyncIndividual(td,i);
        
    //printf("\nfemale: %d, male: %d, offspring: %d to ", f->s, m->s, i->s);
    //printf("offspring ID: %llu, sex: %d, %d of %d\n",i->id,i->s,j+1,rn);
    //if( *(i->f0) == 4 && *(i->f1) == 4 )
    //if(Gen > 135)
    //  printf("p-DRIVE: o(%llu, %llu), m(%llu, %llu), f(%llu, %llu)\n", *(i->f0), *(i->f1), *(m->f0), *(m->f1), *(f->f0), *(f->f1));
      
    //printf("****offspring disperses\n");
    Disperse(td, i, p);
    
    /*
    // sanity print
    if( (*(i->f0) == 1 || *(i->f1) == 1) && (*(f->f0) != 1 && *(f->f1) != 1 && *(m->f0) != 1 && *(m->f1) != 1) ){
      //printf("DRIVE: o(%llu, %llu), m(%llu, %llu), f(%llu, %llu)\n", *(i->m0), *(i->m1), *(m->m0), *(m->m1), *(f->m0), *(f->m1));
      printf("PROLACTIN: o(%llu, %llu), m(%llu, %llu), f(%llu, %llu)\n", *(i->f0), *(i->f1), *(m->f0), *(m->f1), *(f->f0), *(f->f1));
    }
    */
    
    //printf("%d of %d--------\n",j+1,rn);
  }
  //printf("-------------%d\n",rn);
  SinglePaternity++;
  //printf("SinglePaternity: %d--------\n",SinglePaternity);
}


// mates m and f, and there is fertility selection
static void PolyMate(pThreadData td, pIndividual f, pIndividual m_one, pIndividual m_two, pPatch p)
{
  ArbiterProgress();
  pIndividual i;
  pIndividual m;  // this is the siring male
  int         j, rn, paternityindex;
  double      mknockout, fknockout, nc = 0.;
  u64b        l;
  
  //printf("PolyMating attempt btw f: %d(%d), m1: %d(%d), and m2: %d(%d)\n", f->id, f->s, m_one->id, m_one->s, m_two->id, m_two->s);
  //printf("\nPolyMating attempt btw f: %d (sex: %d) with dr: %lf, p: %lf; m1: and %d (sex: %d) with dr: %lf, p: %lf; m2: and %d (sex: %d) with dr: %lf, p: %lf\n", f->id, f->s, f->m, f->f, m_one->id, m_one->s, m_one->m, m_one->f, m_two->id, m_two->s, m_two->m, m_two->f);

  nc = _b*w(td,f,p);
  //nc = w(td,f,p);
  rn = Poisson(&R, nc);
  //printf("no. of offspring seed (w) for Poisson: %f - no. of offspring: %d (%d, %d)\n", nc, rn, p->lat, p->lon);
  OffspringPopSize += rn;
  paternityindex = 0;
  
    
  for(j=0; j<((int)(rn)); j++) {
    //printf("offsp. %d of %d\n", j+1, rn);
    // Get an individual to use as a child
    i = GetNextIndividual(td);
    i->age = 0;
    
    // pick the father here amon n=2 males and initialize to m and continue as before
    int siringmaleinteger = 3;
    siringmaleinteger = rnd(&R, 2);
    //printf("siringmaleinteger: %d\n", siringmaleinteger);
    if( siringmaleinteger == 0 ){
      m = m_one;
    }
    if( siringmaleinteger == 1 ){
      m = m_two;
    }
    
    if( Gen < 174 && ((*(f->f0) == 2 && *(f->f1) == 1) || (*(f->f0) == 1 && *(f->f1) == 2)) ) MatingMotherwithDrive++;
    if( Gen < 174 && ((*(m->f0) == 2 && *(m->f1) == 1) || (*(m->f0) == 1 && *(m->f1) == 2)) ) MatingFatherwithDrive++;
    
    // use this index to calculate multiple paternity if all the father are same, number of offspring should
    // equal to siringmale integer if sired by first male; or should equal to 2Xsiring male integer
    // if sired by the latter. if not equal, then multiple paternity
    paternityindex += siringmaleinteger;
    //printf("paternityindex: %d\n", paternityindex);
    
    // Recombine and mutate
    FastRecombine(td, i, m, f);

    double fgerm = 0.;
    double mgerm = 0.;

    // ***************** drive
    // the following happens in both the approaches (one and two)
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // in _HOMING: the drive is inserted in fertility locus and is happlosufficient
    // in _HOMER:   the drive is inserted in VIABILITY locus and is happloINsufficient
    // in _HOMERKO: the drive is inserted in viability locus and is happloINsufficient
    // in _HOMERKO there is additional knockout in females (or both sexes) of a distal gene
    // if it has the drive (the allele on the other chr should be wildtype (no homing if functional resistant because already cut))
    // has to be heterozygous otherwise regular transmission
    if( (*(f->f0) == 2 && *(f->f1) == 1) || (*(f->f0) == 1 && *(f->f1) == 2) ){
      fgerm = GetFemaleGermLine();
      //if( Gen == 122 ) { printf("female germ: %llf (%llu, %llu)\n", fgerm, *(f->f0), *(f->f1)); }
      // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
      // FEMALES F1
      if( ((int)(fgerm)) == 0 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 2, 3);
      }
      if( ((int)(fgerm)) == 1 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 4, 3);
      }
      if( ((int)(fgerm)) == 2 ){
        // never this in females if pC = 1
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 1, 3);
      }
      if( ((int)(fgerm)) == 3 ){
        // never this in females if pL = 1
        for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 3, 3);
      }
    }
     
    // the following happens in both the approaches (one and two)
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // in _HOMING: the drive is inserted in fertility locus and is happlosufficient
    // in _HOMER:   the drive is inserted in VIABILITY locus and is happloINsufficient
    // in _HOMERKO: the drive is inserted in viability locus and is happloINsufficient
    // in _HOMERKO there is additional knockout in females (or both sexes) of a distal gene
    // if it has the drive (the allele on the other chr should be wildtype (no homing if functional resistant because already cut))
    // has to be heterozygous otherwise regular transmission
    if( (*(m->f0) == 2 && *(m->f1) == 1) || (*(m->f0) == 1 && *(m->f1) == 2) ){
      mgerm = GetMaleGermLine();
      //if( Gen == 122 ) { printf("male germ: %llf (%llu, %llu)\n", mgerm, *(m->f0), *(m->f1)); }
      // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
      // MALES F0
      if( ((int)(mgerm)) == 0 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 2, 3);
      }
      if( ((int)(mgerm)) == 1 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 4, 3);
      }
      if( ((int)(mgerm)) == 2 ){
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 1, 3);
      }
      if( ((int)(mgerm)) == 3 ){
        // never this in males if pL = 1
        for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 3, 3);
      }
    }

#if _HOMER_HS
#if _EXPOSURE
    // the drive allele has changed from 2 to 5, which has lost the rescue, but still "homing"
    if( gen_exp > 0 && Gen >= gen_exp ){
      // get female germline
      if( (*(f->f0) == 5 && *(f->f1) == 1) || (*(f->f0) == 1 && *(f->f1) == 5) ){
        fgerm = GetFemaleGermLine();
        //if( Gen == 122 ) { printf("female germ: %llf (%llu, %llu)\n", fgerm, *(f->f0), *(f->f1)); }
        // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant, 5: drive withour rescue after exposure
        // FEMALES f1
        if( ((int)(fgerm)) == 0 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 5, 3);
        }
        if( ((int)(fgerm)) == 1 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 4, 3);
        }
        if( ((int)(fgerm)) == 2 ){
          // never this in females if pC = 1
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 1, 3);
        }
        if( ((int)(fgerm)) == 3 ){
          // never this in females if pL = 1
          for(l=0; l<_Lf; l++) SetLocus(i->f1, l, 3, 3);
        }
      }
     
      // the following happens in the approaches
      // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
      // in _HOMING:  the drive is inserted in fertility/viability locus and is happlosufficient
      // in _HOMER:   the drive is inserted in VIABILITY locus and is happloINsufficient
      // in _HOMERKO: the drive is inserted in VIABILITY locus and is happloINsufficient
      // in _HOMERKO there is additional knockout in females (or both sexes) of a distal gene
      // if it has the drive (the allele on the other chr should be wildtype (no homing if functional resistant because already cut))
      // has to be heterozygous otherwise regular transmission
      // the following "homing" happens in all three approaches
      if( (*(m->f0) == 5 && *(m->f1) == 1) || (*(m->f0) == 1 && *(m->f1) == 5) ){
        mgerm = GetMaleGermLine();
        //if( Gen == 122 ) { printf("male germ: %llf (%llu, %llu)\n", mgerm, *(m->f0), *(m->f1)); }
        // prolactin alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant, 5: drive withour rescue after exposure
        // MALES f0
        if( ((int)(mgerm)) == 0 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 5, 3);
        }
        if( ((int)(mgerm)) == 1 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 4, 3);
        }
        if( ((int)(mgerm)) == 2 ){
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 1, 3);
        }
        if( ((int)(mgerm)) == 3 ){
          // never this in males if pL = 1
          for(l=0; l<_Lf; l++) SetLocus(i->f0, l, 3, 3);
        }
      }
    }
#endif
#endif

      
    // ***************** knockout in distal loci
    // m locus alleles are 0: nonfunctional, 1: wildtype (_HOMERKO)
    // after the offspring is produced, knockout the genes
    // knockout happens in females in A and C; both sexes in B and D (so also include B and D in the first '#if' conditon below)
    // target for knockout is female fertility in A and B; and male fertility in C and D
    // if both A and C (or both B and D) are both on 'viability' in both sexes (Main())
#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C
    // the fertility knockout is in the germline only (in ALL approach2's!!!)
    // if the mother has the drive, then the offspring will inherit '0' always!
    // ***insert probability of knockout
    // no resistance
    // check if the mother (1st cond) has the drive (2nd cond),
    // THEN knockout fertility (could also check IF it has prolactin, or just knockout anyways)
    // check if the female has the drive in 'f' loci
    if ( (f->s == 1) && ( *(f->f0) == 2 || *(f->f1) == 2 ) ){
      //knockout the 'm' locus in females (assume it inherits 'm1')
      fknockout = 0.;
      fknockout = U01(&R);
      if( fknockout <= _FEMALEKNOCKOUT ){
        for(l=0; l<_Lm; l++) SetLocus(i->m1, l, 0, 1);
      }
    }

    // the fertility knockout is in the male germline as well in B and D in addition to female knockout in A and C
    // if the father has the drive (in this additional option), then the offspring will inherit '0' always!
    // ***insert probability of knockout
    // no resistance
    // check if the father (1st cond) has the drive (2nd cond),
    // THEN knockout fertility (could also check IF it has prolactin, or just knockout anyways)
    // check if the male has the drive in 'f' loci
    if ( (m->s == 0) && ( *(m->f0) == 2 || *(m->f1) == 2 ) ){
      mknockout = 0.;
      mknockout = U01(&R);
      if( mknockout <= _MALEKNOCKOUT ){
        //knockout the 'm' locus in males (assume it inherits 'm0')
        for(l=0; l<_Lm; l++) SetLocus(i->m0, l, 0, 1);
      }
    }
#endif
        
    // ***************** viability check
#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C || _HOMER
    // haploinsufficient: check for viability of offspring in 'f' locus
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // one copy of '4' will make it inviable (allele '2' with drive also has 'rescue')
    if ( *(i->f0) == 4 || *(i->f1) == 4  ){
      // this should work
      continue;
    }
#endif

#if _HOMING_C
    // haploinsufficient: check for viability of offspring in 'f' locus
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    // one copy of '4' will make it inviable (allele '2' with drive also has 'rescue')
      if ( (*(i->f0) == 2 && *(i->f1) == 2) ||
           (*(i->f0) == 2 && *(i->f1) == 4) ||
           (*(i->f0) == 4 && *(i->f1) == 2) ||
           (*(i->f0) == 4 && *(i->f1) == 4) ){
        // this should work
      continue;
    }
#endif

#if _HOMERKO_C
    // // haplosufficient: check for viability of offspring in 'm' locus in C
    // alleles are 1: wildtype, 0: nonfunctional
    if ( *(i->m0) == 0 && *(i->m1) == 0  ){
      // this should work
      continue;
    }
#endif

    // Sync, Handle dispersal and install child
    SyncIndividual(td,i);
    
    //printf("\nfemale: %d, male: %d, offspring: %d to ", f->s, m->s, i->s);
    //printf("offspring ID: %llu, sex: %d, %d of %d\n",i->id,i->s,j+1,rn);
    //if( *(i->f0) != 1 || *(i->f1) != 1 )
      //printf("p-DRIVE: o(%llu, %llu), m(%llu, %llu), f(%llu, %llu)\n", *(i->f0), *(i->f1), *(m->f0), *(m->f1), *(f->f0), *(f->f1));
    //if( *(i->f0) == 4 && *(i->f1) == 4 )
    //if(Gen > 135)
    //  printf("pm-p-DRIVE: o(%llu, %llu), m(%llu, %llu), f(%llu, %llu)\n", *(i->f0), *(i->f1), *(m->f0), *(m->f1), *(f->f0), *(f->f1));
      
    //printf("****offspring disperses in polymate()\n");
    Disperse(td, i, p);
    
    /*
    //sanity print
    if( (*(i->f0) == 1 || *(i->f1) == 1) && (*(f->f0) != 1 && *(f->f1) != 1 && *(m->f0) != 1 && *(m->f1) != 1) ){
      //printf("DRIVE: o(%llu, %llu), m(%llu, %llu), f(%llu, %llu)\n", *(i->m0), *(i->m1), *(m->m0), *(m->m1), *(f->m0), *(f->m1));
      printf("poly-PROLACTIN: o(%llu, %llu), m(%llu, %llu), f(%llu, %llu)\n", *(i->f0), *(i->f1), *(m->f0), *(m->f1), *(f->f0), *(f->f1));
    }
    */
    //printf("%d of %d--------\n",j+1,rn);
  }
  //printf("-------------%d\n",rn);
  if( paternityindex != 0 && paternityindex != rn ){
    MultiplePaternity++;
  }
  //if( Gen >= 14 && MultiplePaternity > 1 ){
    //printf("MultiplePaternity: %d--------\n",MultiplePaternity);
  //}
}


// survival of an individual
static void Survive(pThreadData td, pPatch p, pIndividual m)
{
  ArbiterProgress();
  pIndividual i;
  u64b        l;
  
  // copy all details
  i = GetNextIndividual(td);
  for(l=0; l<T_LENGTH; l++){
    i->d[l] = m->d[l];
    //printf("new: %d, old: %d\n",i->d[l], m->d[l]);
  }
  i->s = m->s;
  i->age = (m->age+1); //increment age, here since copying to Next
  SyncIndividual(td,i);
  //printf("**copying ind - old ID: %d, new ID: %d, sex: %d, age: %d\n", m->id, i->id, i->s, i->age);

  //printf("****adult disperses\n");
  Disperse(td, i, p);
}

/*
#if _EXPOSURE
// survival of an individual
static void Expose(pThreadData td, pPatch p, pIndividual m)
{
  ArbiterProgress();
  pIndividual i;
  u64b        c,l;
  
  // copy all details
  i = GetNextIndividual(td);
  for(l=0; l<T_LENGTH; l++){
    i->d[l] = m->d[l];
    //printf("%d before exposed\n",i->d[l]);
  }
  // f loci: cut out the rescue in Homer_HS so that drive carring individuals become haploinsufficient when 'exposed'
  // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant, 5: drive without rescue
  if( *(m->f0) == 2 ){
    for(c=0; c<_Lf; c++)
        SetLocus(i->f0, c, 5, 3);
  }
  
  if( *(m->f1) == 2 ){
    for(c=0; c<_Lf; c++)
        SetLocus(i->f1, c, 5, 3);
  }

  i->s = m->s;
  i->age = (m->age+1); //increment age, here since copying to Next
  SyncIndividual(td,i);
  //sanity print
  //for(l=0; l<T_LENGTH; l++){
  //  printf("%d after exposed\n",i->d[l]);
  //}
  //printf("**copying ind in exposed - old ID: %d, new ID: %d, sex: %d, age: %d\n", m->id, i->id, i->s, i->age);

  //printf("****adult disperses\n");
  Disperse(td, i, p);
}
#endif
*/

#if _INOC
// mates m and f, and there is fertility selection
static void Inoculate(pThreadData td, pPatch p)
{
  ArbiterProgress();
  pIndividual i;
  int         j;
  u64b        l,c;
  
  for( j=0; j<_INOC_INDIVIDUALS; j++ ) {
    i = GetNextIndividual(td);
    i->s = 0;
    for(l=0; l<T_LENGTH; l++) i->d[l] = 0;
    

    // f loci: drive caryying individuals (either inserted in fertility gene (approachone) or on viability gene (approachtwo))
    // alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
    for(c=0; c<_Lf; c++)     SetLocus(i->f0, c, 2, 3);
    for(c=0; c<_Lf; c++)     SetLocus(i->f1, c, 1, 3);

#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C
    // m loci: distal fertility (or sterility, or, viability) target gene
    // alleles are 0: nonfunctional, 1: drive
    for(c=0; c<_Lm; c++)     SetLocus(i->m0, c, 0, 1);
    for(c=0; c<_Lm; c++)     SetLocus(i->m1, c, 1, 1);
#endif

#if _MALEINOCULATION
    // here _Ln is used as a sex chromosome; heteroz (0,1) are males; homoz. (1,1) are females.
    // gene-drive could be male only, or two sex
    // (e.g. t-allele and x-shredder are both male drives)
    // male only drive
    for(c=0; c<_Ln; c++)     SetLocus(i->z0, c, 0, 4);
    for(c=0; c<_Ln; c++)     SetLocus(i->z1, c, 1, 4);
    i->s = 0;
#endif

#if _FEMALEINOCULATION
    // here _Ln is used as a sex chromosome; heteroz (0,1) are males; homoz. (1,1) are females.
    // gene-drive could be male only, or two sex
    // (e.g. t-allele and x-shredder are both male drives)
    // male only drive
    for(c=0; c<_Ln; c++)     SetLocus(i->z0, c, 1, 4);
    for(c=0; c<_Ln; c++)     SetLocus(i->z1, c, 1, 4);
    i->s = 1;
#endif
 
    i->age = 0;

    SyncIndividual(td,i);
    //printf("sex: %d, f:%lf (%d/%d)\n",i->s, i->f, *i->f0, *i->f1 );
    
    //int n;
    //printf("----inoculating with drive: %llu/%llu and prolactin: %llu/%llu (sex: %d)\n",*(i->m0),*(i->m1),*(i->f0),*(i->f1), i->s);
    //printf("----inoculating with prolactin embedded drive: %llu/%llu (sex: %d)\n",*(i->f0),*(i->f1), i->s);
    //printf("**ind - no: %d, ID: %d, sex: %d\n", j, i->id, i->s);
    //printf("**inoc - no: %d, sex: %d, f: %lf\n", j, i->s, i->f);
    //printf("**inoc - no: %d\n", j);
    //for(n=0; n<_Lf; n++) printf("%llu ",GetLocus(i->f0,n,4));
    //printf("\n");
    //for(n=0; n<_Lf; n++) printf("%llu ",GetLocus(i->f1,n,4));
    //printf("\n");
    //printf("*&*&*&*&*&*&*&*&*&*&\n");
    //for(l=0; l<T_LENGTH; l++) printf("%d ", i[j].d[l]);
    //printf("\n");

    //printf("in inoculate\n");
    Disperse(td, i, p);
  }
}
#endif


// Does all the real work involved in the simulation.
void* THREAD(void *arg)
{

  ThreadData    td;
  u64b          h,i,j,k,l;
  distribution *fdist;
  distribution *pdist;
  pIndividual   secondmatingmale;

  int           g,b;
  pIndividual   matingfemale;
  pIndividual   matingmale;
  pIndividual   female;
  pPatch        matingpatch;
  
  // Initialize this thread's portion of the simulation
  td.id = *((int*)arg);
  Threads[td.id] = &td;
  InitThread(&td);
  
  fdist = allocdist(_MAX_INDIVIDUALS*_DRANGE*_DRANGE);
  pdist = allocdist(_MAX_INDIVIDUALS*_DRANGE*_DRANGE);

  while(SignalDone(&td)){

#if _DISTANCE
    for( i=0; i<12; i++){
      for( j=0; j<((int)(_DRANGE+1)); j++){
        distancedata[i][j] = 0;
        //printf("%d ", distancedata[i][j]);
      }
    }
#endif
    
    //----------------
#if _EXPOSURE
    u64b        c,l;
    for(i=0; i<WIDTH; i++){
      //printf("\n");
      for(j=0; j<HEIGHT; j++){
        for(h=0, l=td.Current[i][j].ni; h<l; h++){          // loops through individuals in a patch
          int exposureflag = 0;
            double randomexp = 0.;                          // exposure probability
            randomexp = U01(&R);
            if( gen_exp > 0 && Gen == gen_exp && randomexp <= _ExposureProbability ){
            if( *(td.Current[i][j].i[h]->f0) == 2 || *(td.Current[i][j].i[h]->f1) == 2 ){
              exposureflag = 1;
            }
          }
              
          if( exposureflag == 1 ){
            //printf("before exp: %d, %d, age: %d, ([%d, %d] %d/%d)\n", *(td.Current[i][j].i[h]->f0), *(td.Current[i][j].i[h]->f1), td.Current[i][j].i[h]->age, i, j, h+1, l);
            if( *(td.Current[i][j].i[h]->f0) == 2 ){
              for(c=0; c<_Lf; c++) SetLocus(td.Current[i][j].i[h]->f0, c, 5, 3);
            }
              
            if( *(td.Current[i][j].i[h]->f1) == 2 ){
              for(c=0; c<_Lf; c++) SetLocus(td.Current[i][j].i[h]->f1, c, 5, 3);
            }

              //Expose(&td,&td.Next[i][j],td.Current[i][j].i[h]);
            //printf("after exp: %d, %d, age: %d, ([%d, %d] %d/%d)\n", *(td.Current[i][j].i[h]->f0), *(td.Current[i][j].i[h]->f1), td.Current[i][j].i[h]->age, i, j, h+1, l);
          }
        }
      }
    }
#endif
    //----------------
    
      
    OffspringPopSize = 0;
    Mortality = 0;
    SinglePaternity = 0;
    MultiplePaternity = 0;

    //printf("in thread.c %lf\n",_LOFPROB);
    for(i=0; i<WIDTH; i++){
      //printf("\n");
      for(j=0; j<HEIGHT; j++){
        int numberoffemales = 0;
        int numberofmales   = 0;
        int mating          = 0;
        int success         = 0;
        int maleinteger     = 0;
        int nothers         = 0;
        double maleprobsum  = 0.;
            
#if _DISTANCEDENSITYDISPERSAL
        td.Current[i][j].density = ( ((double)(NeDispersal(&td,&td.Current[i][j]))) / ((double)(td.Current[i][j].ddne * _K0)) );
        //printf("before-----(%d, %d) N: %d, D: %lf \n", i,j,td.Current[i][j].ni, td.Current[i][j].density);
#endif
        //if (td.Current[i][j].ni > 0) printf("-----(%d, %d) N: %d \n", i,j,td.Current[i][j].ni);
        
        for(h=0, l=td.Current[i][j].ni; h<l; h++){          // loops through individuals in a patch
          //printf("%d of %d, sex: %d\n",h,l, td.Current[i][j].i[h]->s);
          //printf("%d of %d, sex: %d, f:%lf (%llu/%llu)\n",h,l, td.Current[i][j].i[h]->s, td.Current[i][j].i[h]->f, *td.Current[i][j].i[h]->f0, *td.Current[i][j].i[h]->f1 );
          if( td.Current[i][j].i[h]->s == MALE ) continue;  // Skip any males we encounter
          female = td.Current[i][j].i[h];                   // the individual is a female, find a mate:
          
#if _HOMING_A || _HOMING_C
          // skip infertile females; haplosufficient, both copies should be missing
          // printf("f0: %llu, f1: %llu\n", *(female->f0), *(female->f1));
          if( *(female->f0) == 2 && *(female->f1) == 2 ) continue;
          if( *(female->f0) == 2 && *(female->f1) == 4 ) continue;
          if( *(female->f0) == 4 && *(female->f1) == 2 ) continue;
          if( *(female->f0) == 4 && *(female->f1) == 4 ) continue;
#endif

#if _HOMERKO_A || _HOMERKO_C
          // m locus alleles are 0: nonfunctional, 1: drive (_HOMERKOA and _HOMERKOC only)
          // females are infertile if both copies are missing
          if( *(female->m0) == 0 && *(female->m1) == 0 ) continue;
#endif

#if _HOMER_HS
          // females are infertile if both copies are missing
          if( *(female->f0) == 4 && *(female->f1) == 4 ) continue;
#if _EXPOSURE
          if( gen_exp > 0 && Gen >= gen_exp ){
            if( ( *(female->f0) == 5 && *(female->f1) == 5  ) ||
                ( *(female->f0) == 4 && *(female->f1) == 5  ) ||
                ( *(female->f0) == 5 && *(female->f1) == 4  ) ){
              continue;
            }
          }
#endif
#endif

            
            
          //printf("----choosing: %llu/%llu\n",*(female->f0),*(female->f1));
          matingfemale = td.Current[i][j].i[h];             // the female is fertile, find a mate:
          success = 0;
          numberoffemales++;
          //printf("(%d, %d) mating female (no: %d, %d) searching...... \n", i,j, matingfemale->id, matingfemale);
          //printf("\nfemale searching...... \n");
          
          ArbiterProgress();                                // this is for Aaron's arbiter
          //--------------------------------------------------------------------------------------
          //--------------------------------------------------------------------------------------
          //--------------------------------------------------------------------------------------

#if _PANMICTIC
          // start a random array and go through it systematically
          // this will be critical when pop sizes are low
          int randomarray[(_WIDTH*_HEIGHT)],v,temp,randomIndex;
          for( v = 0; v < (_WIDTH*_HEIGHT); v++ ) { randomarray[v] = v; }
                    
          for( v = 0; v < (_WIDTH*_HEIGHT); v++ ) {
            temp = randomarray[v];
            randomIndex = rnd(&R, (_WIDTH*_HEIGHT));
            randomarray[v] = randomarray[randomIndex];
            randomarray[randomIndex] = temp;
          }

          // find a patch randomly
          int rlat = 0;
          int rlon = 0;
          int mg;
          // loop through as many times as there are patches until a male is found
          for(mg=0; mg < (_WIDTH*_HEIGHT); mg++ ){
            rlat = randomarray[mg]/(_WIDTH);
            rlon = randomarray[mg]%(_WIDTH);
                      
            if( !td.Current[rlat][rlon].ni) continue;         // if no ind. here in the random patch skip it...
            matingpatch = &td.Current[rlat][rlon];
                      
            // there is someone in the patch, loop through individuals
            int gn;
            numberofmales = 0;
            maleprobsum = 0.0;
            for(gn=0; gn<matingpatch->ni; gn++ ){
              if( matingpatch->i[gn]->s == MALE ){
                numberofmales++;
                maleprobsum += ( fdist->p[gn] = 1. );
              } else {
                fdist->p[gn] = 0.;
              }
            }

            if( !numberofmales || !maleprobsum ) continue;

            // if there are males in this patch, break the search,then mate...
            if( maleprobsum != 0. ){
              maleinteger = 0;
              fdist->n = matingpatch->ni;
              initdist(fdist,maleprobsum);
              maleinteger = drand(&R, fdist);
              matingmale = matingpatch->i[maleinteger];
              success = 1;    // initialize the success to 1, so that we don't force mating with maleintr = 0;
              break;
            }
          } // closes the gn loop that goes through all the patches systematically

#if _HOMING_B || _HOMING_C
          if( (*(matingmale->f0) == 2 && *(matingmale->f1) == 2) ||
              (*(matingmale->f0) == 2 && *(matingmale->f1) == 4) ||
              (*(matingmale->f0) == 4 && *(matingmale->f1) == 2) ||
              (*(matingmale->f0) == 4 && *(matingmale->f1) == 4) ){
            success = 0;
          }
#endif
#if _EXP_MALE_STERILITY
          if( *(matingmale->f0) == 2 || *(matingmale->f1) == 2 ) {
            success = 0;
          }
#endif
#if _HOMERKO_B || _HOMERKO_C
          if( *(matingmale->m0) == 0 && *(matingmale->m1) == 0 ) {
            success = 0;
          }
#endif

          if( success != 0 ){
            //printf("panmictic mating\n");
            Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
            mating++;
          }
          // the global search should be over by now...
          //--------------------------------------------------------------------------------------
          //--------------------------------------------------------------------------------------
          //--------------------------------------------------------------------------------------

#else
          // SEARCH 1
          // searching at the 'home' patch
          int n;
          numberofmales = 0;
          maleprobsum   = 0.0;
          for(n=0; n<td.Current[i][j].ni; n++ ){
            if( td.Current[i][j].i[n]->s == MALE ){
              numberofmales++;
              //printf("number males: %d(%d, %d)\n", numberofmales, i, j);
              maleprobsum += ( fdist->p[n] = 1. );
              pdist->p[n] = 1.;
            } else {
              fdist->p[n] = 0.0;
              pdist->p[n] = 0.0;
            }
            //printf("fdist: %lf\n", fdist->p[n]);
          }
          //if (numberofmales == 0 ) printf("no males at home\n");
          //printf("maleprobsum: %lf\n", maleprobsum);
          
          if( maleprobsum != 0. ){
            maleinteger = 0;
            fdist->n = td.Current[i][j].ni;
            initdist(fdist,maleprobsum);
            maleinteger = drand(&R, fdist);
            matingmale = td.Current[i][j].i[maleinteger];
            success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
          }
 
          double polymatingprob = 0.;
          int secondmaleinteger = 0;
          int polymating        = 0;
          int polysuccess       = 0;
          
          polymatingprob = U01(&R);
          //printf("polymatingprob: %lf\n", polymatingprob);
          // if more than 1 male, the first male is chosen with success and prob is satisfied, then:
          if( numberofmales > 1 && success == 1 && polymatingprob <= _POLYANDRYPROB) {
            //printf("[%d] %lf", maleinteger, pdist->p[maleinteger]);
            pdist->p[maleinteger] = 0.0;
            //printf("after [%d] %lf", maleinteger, pdist->p[maleinteger]);
            pdist->n = td.Current[i][j].ni;
            initdist(pdist,(maleprobsum-1));
            secondmaleinteger = drand(&R, pdist);
            //printf("numberofmales: %d, maleprobsum: %lf, polymatingprob: %lf, maleinteger: %d, secondmaleinteger: %d\n", numberofmales, maleprobsum, polymatingprob, maleinteger, secondmaleinteger);
            if( (secondmaleinteger - maleinteger) == 0){
              printf("warning! picking the same male twice at (%d, %d): %d - %d\n", i, j, maleinteger, secondmaleinteger);
            }
            secondmatingmale = td.Current[i][j].i[secondmaleinteger];
            
            polysuccess = 1;
#if _HOMING_B || _HOMING_C
            if( (*(secondmatingmale->f0) == 2 && *(secondmatingmale->f1) == 2) ||
                (*(secondmatingmale->f0) == 2 && *(secondmatingmale->f1) == 4) ||
                (*(secondmatingmale->f0) == 4 && *(secondmatingmale->f1) == 2) ||
                (*(secondmatingmale->f0) == 4 && *(secondmatingmale->f1) == 4) ){
              polysuccess = 0;
            }
#endif
#if _EXP_MALE_STERILITY
            if( *(secondmatingmale->f0) == 2 || *(secondmatingmale->f1) == 2 ) {
              polysuccess = 0;
            }
#endif
#if _HOMERKO_B || _HOMERKO_C
            // this is a bit not so straightforward
            // male fertility is affected in C and D only
            // BUT 'viability' is implied in 'm' locus if both ("A&C" are on simultaneously; or "B&D" are on simultaneously)
            // target is 'm' loci; males are infertile if both copies are missing
            // m locus alleles are 0: nonfunctional, 1: drive
            if( (*(secondmatingmale->m0) == 0 && *(secondmatingmale->m1) == 0) ){
              polysuccess = 0;
            }
#endif
          }

          //--------------------------------------------------------------

          // originally we skip these males (see below), but if polyandry, and the second male she pick is ok
          // no need to continue, but Mate() with the second male only
          int newpolytosingleflag = 0;
#if _HOMING_B || _HOMING_C
          if ( polysuccess != 0 ){
            // skip infertile males; haplosufficient, both copies should be missing
            // printf("f0: %llu, f1: %llu\n", *(female->f0), *(female->f1));
            if( *(matingmale->f0) == 2 && *(matingmale->f1) == 2 ) { newpolytosingleflag = 1; }
            if( *(matingmale->f0) == 2 && *(matingmale->f1) == 4 ) { newpolytosingleflag = 1; }
            if( *(matingmale->f0) == 4 && *(matingmale->f1) == 2 ) { newpolytosingleflag = 1; }
            if( *(matingmale->f0) == 4 && *(matingmale->f1) == 4 ) { newpolytosingleflag = 1; }
            if( newpolytosingleflag == 1 ) {
              matingmale = secondmatingmale;
              newpolytosingleflag = 0;
            }
          }
#endif
#if _EXP_MALE_STERILITY
          if ( polysuccess != 0 ){
            // skip infertile males; haplosufficient, both copies should be missing
            // printf("f0: %llu, f1: %llu\n", *(female->f0), *(female->f1));
            if( *(matingmale->f0) == 2 || *(matingmale->f1) == 2 ) { newpolytosingleflag = 1; }
            if( newpolytosingleflag == 1 ) {
              matingmale = secondmatingmale;
              newpolytosingleflag = 0;
            }
          }
#endif
#if _HOMERKO_B || _HOMERKO_C
          // this is a bit not so straightforward
          // male fertility is affected in C and D only
          // BUT 'viability' is implied in 'm' locus if both ("A&C" are on simultaneously; or "B&D" are on simultaneously)
          // target is 'm' loci; males are infertile if both copies are missing
          // m locus alleles are 0: nonfunctional, 1: drive
          if ( polysuccess != 0 ){
            // skip infertile males; haplosufficient, both copies should be missing
            // printf("f0: %llu, f1: %llu\n", *(female->f0), *(female->f1));
            if( *(matingmale->m0) == 0 && *(matingmale->m1) == 0 ) { newpolytosingleflag = 1; }
            if( newpolytosingleflag == 1 ) {
              matingmale = secondmatingmale;
              newpolytosingleflag = 0;
            }
          }
            
#endif
          //--------------------------------------------------------------
          if( newpolytosingleflag == 1 || polysuccess == 0){
            if( success != 0 ){
#if _HOMING_B || _HOMING_C
              if( (*(matingmale->f0) == 2 && *(matingmale->f1) == 2) ||
                  (*(matingmale->f0) == 2 && *(matingmale->f1) == 4) ||
                  (*(matingmale->f0) == 4 && *(matingmale->f1) == 2) ||
                  (*(matingmale->f0) == 4 && *(matingmale->f1) == 4) ){
                continue;
              }
#endif
#if _EXP_MALE_STERILITY
              if( *(matingmale->f0) == 2 || *(matingmale->f1) == 2 ) {
                continue;
              }
#endif
#if _HOMERKO_B || _HOMERKO_C
              if( *(matingmale->m0) == 0 && *(matingmale->m1) == 0 ) {
                continue;
              }
#endif

              //if( *(matingmale->f0) != 1 || *(matingmale->f1) != 1 )
                //printf("mating at 'home': m(%llu, %llu)\n", *(matingmale->f0), *(matingmale->f1));
                
              //printf("mating at 'home': f: %d (%d), m: %d (%d)\n",matingfemale, matingfemale->s, matingmale, matingmale->s);
              Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
              //printf("...mated at 'home'\n");
              mating++;
              //continue; // (no need i think, so silenced) continue so that she doesn't mate multiple times
            }
          } else if ( newpolytosingleflag == 0 && polysuccess == 1) {
            //if( *(matingmale->f0) != 1 || *(matingmale->f1) != 1 )
              //printf("p mating at 'home': m1(%llu, %llu), m2(%llu, %llu)\n", *(matingmale->f0), *(matingmale->f1), *(secondmatingmale->f0), *(secondmatingmale->f1) );
              
            PolyMate(&td,matingfemale,matingmale,secondmatingmale,&td.Next[i][j]);
            //printf("...poly. mated at 'home'\n");
            polymating++;
            //continue; // (no need i think, so silenced) continue so that she doesn't mate multiple times
          }
          // this shouldn't happen:
          //if( success == 1 ) continue;
          
          // SEARCH 2
          // if no males in 'home' patch, search among the roamers at the 'home' patch
          int     numberofmalesroam   = 0;
          double  maleprobsumroam     = 0.;
          success = 0;
          if( !maleprobsum || !numberofmales ){
            int n;
            nothers = td.Current[i][j].nothers;
            if( !nothers) continue;
            //printf("checking roaming inds. no. of others in the patch (%d, %d) is: %d\n", i, j, nothers);
            for(n=0; n<nothers; n++ ){
              if( td.Current[i][j].Other[n].ind->s == MALE ){
                numberofmalesroam++;
                maleprobsumroam += ( fdist->p[n] = 1. );
              } else {
                fdist->p[n] = 0.;
              }
            }
            //printf("N roaming males: %d, P: %lf\n", numberofmalesroam, maleprobsumroam);
          }
        
          //if (numberofmalesroam == 0 ) printf("no males roaming at home\n");
            
          if( maleprobsumroam != 0. && numberofmalesroam !=0 ){
            maleinteger = 0;
            fdist->n = nothers;
            initdist(fdist,maleprobsumroam);
            maleinteger = drand(&R, fdist);
            matingmale = td.Current[i][j].Other[maleinteger].ind;
            success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
            //printf("male integer: %d, male: %d, sex: %d\n", maleinteger, matingmale, matingmale->s);
          }
          
          // skip t-allele homozygous male if picked, the female loses a mating opportunity
          if( success != 0 ){
#if _HOMING_B || _HOMING_C
            if( (*(matingmale->f0) == 2 && *(matingmale->f1) == 2) ||
                (*(matingmale->f0) == 2 && *(matingmale->f1) == 4) ||
                (*(matingmale->f0) == 4 && *(matingmale->f1) == 2) ||
                (*(matingmale->f0) == 4 && *(matingmale->f1) == 4) ){
              continue;
            }
#endif
#if _EXP_MALE_STERILITY
            if( *(matingmale->f0) == 2 || *(matingmale->f1) == 2 ) {
              continue;
            }
#endif
#if _HOMERKO_B || _HOMERKO_C
            if( *(matingmale->m0) == 0 && *(matingmale->m1) == 0 ) {
              continue;
            }
#endif

            //if( *(matingmale->f0) != 1 || *(matingmale->f1) != 1 )
              //printf("mating at 'home' roamer: m(%llu, %llu)\n", *(matingmale->f0), *(matingmale->f1));
              
            //printf("mating at 'home' roamer f: %d (%d), m: %d (%d)\n",matingfemale, matingfemale->s, matingmale, matingmale->s);
            Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
            //printf(".....mated at 'home' with a roaming male\n");
            mating++;
            //continue; // (no need i think, so silenced)
          }
          
          //***************************************
          // SEARCH 3
          // if no males in 'home' patch (+ no roamers)
          int distmatingsuccess;
          
          if( !numberofmales && !success ){
                        
            //int mr; u64b w; double m;
            int matedist,nPatches,totalPatches,matei,matej,v;
            distmatingsuccess = 0;
            // skips matedist=0, should have already found it if it was there.
            // if no males in 'home' patch, search among the roamers at the 'home' patch
          
            //for( matedist=1; matedist<((_RRANGE-1)/2+1); matedist++ ) {
            for( matedist=1; matedist<((_DRANGE-1)/2+1); matedist++ ) {
              
              if( distmatingsuccess == 1 ) continue;
              //printf("distance: %d\n", matedist);
              // for a given distance initiate an array for coordinates to randomly pick from
              
              totalPatches = (matedist*8);
              //int coordArray[_RRANGE*_RRANGE][2];
              int coordArray[_DRANGE*_DRANGE][2];
              int tempArray[1][2];
              int randomPIndex;
              nPatches = 0;
              
              for( matei=-matedist; matei<(matedist+1); matei++ ) {
                for( matej=-matedist; matej<(matedist+1); matej++ ) {
                  if( (!matei) && (!matej) )                continue;
                  if( !ValidMatePatch(&td,i,j,matei,matej)) continue;
                  coordArray[nPatches][0] = i+matei;
                  coordArray[nPatches][1] = j+matej;
                  nPatches++;
                }
              }
              //printf("%d\n",nPatches);
              //for( v = 0; v < (nPatches); v++ ) printf("%d: %d,%d\n",v, coordArray[v][0],coordArray[v][1]);

              // or randomly rearrange the array and go through systematically
              for( v = 0; v < (nPatches); v++ ) {
                tempArray[0][0] = coordArray[v][0];
                tempArray[0][1] = coordArray[v][1];
                randomPIndex = rnd(&R, nPatches);
                coordArray[v][0] = coordArray[randomPIndex][0];
                coordArray[v][1] = coordArray[randomPIndex][1];
                coordArray[randomPIndex][0] = tempArray[0][0];
                coordArray[randomPIndex][1] = tempArray[0][1];
              }
              //for( v = 0; v < (nPatches); v++ ) printf("%d: %d,%d\n",v, coordArray[v][0],coordArray[v][1]);
            
              
              int mg;
              int rlat = 0, rlon = 0, success = 0, maleintr = 0;
              // loop through as many times as there are patches until a male is found
              
              for(mg=0; mg < (nPatches); mg++ ){
                // random (lat, long)
                rlat = coordArray[mg][0];
                rlon = coordArray[mg][1];
                //printf("mg: %d, random: %d, %d\n", mg, rlat, rlon);

                if( !td.Current[rlat][rlon].ni) continue;         // if no ind. here in the random patch skip it
                matingpatch = &td.Current[rlat][rlon];
                //printf("mating patch's (%d,%d) N: %d\n", rlat, rlon, matingpatch->ni);
        
                // there is someone in the patch, loop through individuals
                int gn;
                numberofmales = 0;
                maleprobsum   = 0.0;
                for(gn=0; gn<matingpatch->ni; gn++ ){
                  if( matingpatch->i[gn]->s == MALE ){
                    numberofmales++;
                    maleprobsum += ( fdist->p[gn] = 1. );
                  } else {
                    fdist->p[gn] = 0.;
                  }
                } // closes the gn loop
                //printf("number of males at incremental s: %d\n", numberofmales);
                //if (numberofmales == 0 ) printf("no males at incremental search\n");
                  
                if( maleprobsum != 0. && distmatingsuccess != 1 ){
                  maleinteger = 0;
                  fdist->n = matingpatch->ni;
                  initdist(fdist,maleprobsum);
                  maleinteger = drand(&R, fdist);
                  matingmale = matingpatch->i[maleinteger];
                  success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
                }
  
                // skip t-allele homozygous male if picked, the female loses a mating opportunity
                if( success != 0 && distmatingsuccess != 1){
#if _HOMING_B || _HOMING_C
                  if( (*(matingmale->f0) == 2 && *(matingmale->f1) == 2) ||
                      (*(matingmale->f0) == 2 && *(matingmale->f1) == 4) ||
                      (*(matingmale->f0) == 4 && *(matingmale->f1) == 2) ||
                      (*(matingmale->f0) == 4 && *(matingmale->f1) == 4) ){
                    continue;
                  }
#endif
#if _EXP_MALE_STERILITY
                  if( *(matingmale->f0) == 2 || *(matingmale->f1) == 2 ) {
                    continue;
                  }
#endif
#if _HOMERKO_B || _HOMERKO_C
                  if( *(matingmale->m0) == 0 && *(matingmale->m1) == 0 ) {
                    continue;
                  }
#endif
                  //if( *(matingmale->f0) != 1 || *(matingmale->f1) != 1 )
                    //printf("mating after i. search: m(%llu, %llu)\n", *(matingmale->f0), *(matingmale->f1));
                    
                  //printf("inc. s. mating: f: %d (id: %d) m: %d (id: %d)\n\n",matingfemale, matingfemale->id, matingmale, matingmale->id);
                  Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
                  //printf(".....mated after inc. search.\n");
                  distmatingsuccess = 1;
                  mating++;
                }
                //*(*(*(*(*(*(*(*(
                // SEARCH 4
                // check the roaming males as well here
                if( !distmatingsuccess ){
                  int     n, nothers;
                  int     numberofmalesroam = 0;
                  double  maleprobsumroam   = 0.;
                  
                  nothers = matingpatch->nothers;
                  if( !nothers ) continue;
                  //printf("checking roaming inds. no. of others in the patch (%d, %d) is: %d\n", rlat, rlon, nothers);
                  for(n=0; n<nothers; n++ ){
                    if( matingpatch->Other[n].ind->s == MALE ){
                      //printf("found at least 1 roaming male\n");
                      numberofmalesroam++;
                      maleprobsumroam += ( fdist->p[n] = 1. );
                    } else {
                      fdist->p[n] = 0.;
                    }
                  }
                
                      
                  //printf("number of males: %d, prob: %lf\n", numberofmalesroam, maleprobsumroam);
                  //if (numberofmalesroam == 0 ) printf("no males roaming at incremental search \n");
                  if( !numberofmalesroam || !maleprobsumroam ) continue; // if no males, go to the next patch
                                                
                  if( maleprobsumroam != 0. ){ // if there are males then mate...
                    maleinteger = 0;
                    fdist->n = nothers;
                    initdist(fdist,maleprobsumroam);
                    maleinteger = drand(&R, fdist);
                    matingmale = matingpatch->Other[maleinteger].ind;
                    success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
                  }
  
                  // skip t-allele homozygous male if picked, the female loses a mating opportunity
                  if( success != 0 ){
#if _HOMING_B || _HOMING_C
                    if( (*(matingmale->f0) == 2 && *(matingmale->f1) == 2) ||
                        (*(matingmale->f0) == 2 && *(matingmale->f1) == 4) ||
                        (*(matingmale->f0) == 4 && *(matingmale->f1) == 2) ||
                        (*(matingmale->f0) == 4 && *(matingmale->f1) == 4) ){
                      continue;
                    }
#endif
#if _EXP_MALE_STERILITY
                    if( *(matingmale->f0) == 2 || *(matingmale->f1) == 2 ) {
                      continue;
                    }
#endif
#if _HOMERKO_B || _HOMERKO_C
                    if( *(matingmale->m0) == 0 && *(matingmale->m1) == 0 ) {
                      continue;
                    }
#endif

                    //if( *(matingmale->f0) != 1 || *(matingmale->f1) != 1 )
                      //printf("mating search romaer 'home': m(%llu, %llu)\n", *(matingmale->f0), *(matingmale->f1));
                      
                    //printf("incr. s. roamer mating: f: %d, (%d) m: %d (%d)\n",matingfemale, matingfemale->s, matingmale, matingmale->s);
                    Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
                    //printf(".....mated afters search with a roaming male\n");
                    distmatingsuccess = 1;
                    mating++;
                    continue; // break so that she doesn't mate multiple times within mating 'rings'
                  }
                }
                //*(*(*(*(*(*(*(*(
              } // closes mg loop
              
            } // closes mate dist loop
          } // closes if no males in 'home' patch (+roamers)
#endif
        } // closes the female search loop (going through all N in patch)
        

        for(h=0, l=td.Current[i][j].ni; h<l; h++){
          //printf("\n*checking survival-");
          //printf("(%d, %d) %d of %d\n", i,j,h+1,l);
          if( !l ) continue;                                          // skip empty patches
          if( td.Current[i][j].i[h]->age > _MAXAGE ) continue;        // Skip any inds older than maximum age
        
          // age dependent survival, adults (!=0) have a fixed survival probability
          // this is an attempt to make the overlapping generations less overlap
          double randomno = 0.;
          randomno = U01(&R);
          if( randomno < _SURVIVALPROBABILITY ) {
            //printf("survives: p(s): %lf, random: %lf\n",_SURVIVALPROBABILITY, randomno);
/*
#if _EXPOSURE
            // if the individual is surviving, check if genotype has to be reset for removal of 'rescue'
            // any individuals any '2' should be '5'
            int exposureflag = 0;
            if( gen_exp > 0 && Gen >= gen_exp ){
              if( *(td.Current[i][j].i[h]->f0) == 2 || *(td.Current[i][j].i[h]->f1) == 2 ){
                exposureflag = 1;
              }
            }
            
            if( exposureflag == 1 ){
              Expose(&td,&td.Next[i][j],td.Current[i][j].i[h]);
              //if( (*(td.Current[i][j].i[h]->f0) == 2 || *(td.Current[i][j].i[h]->f1) == 2 ) && gen_exp > 0 && Gen >= gen_exp ){ printf("exp: %d, %d, age: %d, ([%d, %d] %d/%d)\n", *(td.Current[i][j].i[h]->f0), *(td.Current[i][j].i[h]->f1), td.Current[i][j].i[h]->age, i, j, h+1, l);}
            } else {
              Survive(&td,&td.Next[i][j],td.Current[i][j].i[h]);
                if( (*(td.Current[i][j].i[h]->f0) == 2 || *(td.Current[i][j].i[h]->f1) == 2 ) && gen_exp > 0 && Gen >= gen_exp ){ printf("surv: %d, %d, age: %d, ([%d, %d] %d/%d)\n", *(td.Current[i][j].i[h]->f0), *(td.Current[i][j].i[h]->f1), td.Current[i][j].i[h]->age, i, j, h+1, l);}
            }
          } else {
            //printf("dies: p(s): %lf, random: %lf\n",_SURVIVALPROBABILITY, randomno);
              //if( (*(td.Current[i][j].i[h]->f0) == 2 || *(td.Current[i][j].i[h]->f1) == 2 ) && gen_exp > 0 && Gen >= gen_exp ){ printf("dies: %d, %d, age: %d, ([%d, %d] %d/%d)\n", *(td.Current[i][j].i[h]->f0), *(td.Current[i][j].i[h]->f1), td.Current[i][j].i[h]->age, i, j, h+1, l);}
            Mortality++;
          }
#else*/
            Survive(&td,&td.Next[i][j],td.Current[i][j].i[h]);
          } else {
            //printf("dies: p(s): %lf, random: %lf\n",_SURVIVALPROBABILITY, randomno);
            Mortality++;
          }
//#endif
        }
      } // closes j loop (new)
    } //closes i loop (new)
            
    for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
        //printf("(%d, %d) current: %llu, next: %llu\n", i, j, td.Current[i][j].ni, td.Next[i][j].ni);
        td.Current[i][j].ni = 0;
      }
    }

    // Free and swap to generations to make all next gen current.
    ReleaseAllCurrentIndividuals(&td);
    SwapGenerations(&td);
#if _DISTANCE
    LogDistance();
#endif
    // after swapping, next is current, so I can check for survival here
    for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
#if _DISTANCEDENSITYDISPERSAL
        // critical for survival calculations!
        td.Current[i][j].density = ( ((double)(NeDispersal(&td,&td.Current[i][j]))) / ((double)(td.Current[i][j].ddne * _K0)) );
        //printf("after-----(%d, %d) N: %d, D: %lf \n", i,j,td.Current[i][j].ni, td.Current[i][j].density);
#endif
  
#if _INOC
        // this should move to 'after' survival
        for( g = 0; g < _NoInoc; g++ ) {
          if( Gen == InocGen[g] ) {
            
            for( b=0; b < _InocPatchNo; b++ ) {
              if( i == InocCoor[b][0] && j == InocCoor[b][1]) {
                Inoculate(&td,&td.Next[i][j]);
              }
            }
             
            /*
              if( i <_WIDTH && j < _HEIGHT) {
                Inoculate(&td,&td.Next[i][j]);
              }
            */
          }
        }
#endif
      }
    }
    
    for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
        td.Current[i][j].nothers=0;
      }
    }
    
    // calculates who is using which patch for how long
    for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
        for(k=0,l=td.Current[i][j].ni; k < l; k++) {
            TimeSpentRoam(td.Current[i][j].i[k], &td.Current[i][j]);
        }
      }
    }
            
  } // while signal done
  return NULL;
} // void main
