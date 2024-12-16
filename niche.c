#ifndef __NICHE_C
#define __NICHE_C

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

// Fill in an array of thread start functions depending on _CPUs
#if    _CPUS == 1
extern void* Thread_0(void*);
WorkThread WorkThreads[_CPUS] = {Thread_0};
#else
#error !! This code needs to be modified to support _CPUS != 1.
#endif

RndState    R;

//extern void LogHabitat(int clusterno, int *clusterinfo);
//extern void InitThread(pThreadData td);

// Muticies, condition variables, and flags for thread syncronization
volatile    pthread_cond_t  Master,Slaves;
volatile    pthread_mutex_t Mutex;
volatile    int             Go[_CPUS],Done[_CPUS];

// Each thread gets a ThreadData structure which will hold information
// about the thread as well as a copy of some parameters.
// This is the main() thread's ThreadData structure:
ThreadData  Mtd;
pThreadData Mt = &Mtd;

// Array of pointers to the work threads ThreadData strcutures.
// These are malloc'd by the threads themselves.
volatile    pThreadData     Threads[_CPUS];

// Random Seed, and starting space state.  These are filled in by the main
// thread and copied by the work threads.
u64b        Seed;
Space       Start;
int         Map[_WIDTH*_HEIGHT];
// Current generation number and the number of generations to run for
u64b        Gen,NGen;

// Array of all possible niches
Patch       Niches[NICHES];

double      Xi[] = _Xi;
int         MaxPopSize;
//double      _LOFPROB = 1-(1/((double)(_PL)));
int         OffspringPopSize;
int         Mortality;
int         SinglePaternity;
int         MultiplePaternity;
int         MatingMotherwithDrive;
int         MatingFatherwithDrive;

#if _HOMER
int         gen_n;
int         gen_nf;
int         gen_nn;
#endif
#if _EXPOSURE
u64b        gen_exp; //exposure generation when drive frequency exceed _ExposureFrequency
#endif

// log files
#if _SYS
u64b        Sys_logs;
FILE       *Sysf;
#endif
#if _DEME
u64b        Deme_logs;
FILE       *Demef;
#endif
#if _DISTANCE
u64b        Distance_logs;
FILE       *Distancef;
#endif
#if _IND
u64b        Ind_logs;
FILE       *Indf;
#endif
// Benchmark counters
#if _BENCHMARK
u64b        ProcessTime,nProcessed;
#endif

// Code shared with thread.c
#include "shared.h"

static void error( int x, char *s, ... )
{
  va_list args;
  va_start( args, s );
  vfprintf( stderr, s, args );
  if(x)
    exit(x);
}

// A number of functions are placed here which are very similar to 
// those in shared.h.  However, they are only needed in niche.c
// and this will keep the compiler quiet.  You'll find these below:

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


// Returns the number of individuals used in the current
// Note: This count _will_ include wanderers if there are any.
static u64b nCurrentIndividuals()
{
  u64b i,m;

  for(m=i=0; i<_CPUS; i++)
    m += Threads[i]->nICurrent;
  return m;
}


// Syncs the "niche value" of patch p from it's traits
static void SyncPatch(pPatch p)
{
  u64b k;

  for(p->n=k=0; k<_k; k++){
    p->n |= ((u64b)p->l[k])<<k;
  //printf("\n %d---", p->n);
  }
  p->n = log(p->n)/log(2);
  //printf("%d", p->n);
}


// Stop-Watch for temporary timing:
// I could be wrong, but I think all this stopwatch code is only used by the 
// benchmarking code.  (see _BENCHMARK)
struct timeval sw_st,sw_et;

#if 0
// Starts the stopwatch
static void swStart()
{
  gettimeofday(&sw_st,NULL);
}


// Returns the amount of time that has passed since the stopwatch was started.
// So the name of swStop() isn't really the best I suppose.
static double swStop()
{
  gettimeofday(&sw_et,NULL);
  return ((sw_et.tv_sec*((u64b)1000000))+sw_et.tv_usec) -
         ((sw_st.tv_sec*((u64b)1000000))+sw_st.tv_usec);
}
#endif


// Sets m to be the max value in the size l array d, and sets w to be the index of m in d
static void Max(double *d, u64b l, double *m, u64b *w)
{
  int i,c,t[1024];

  // Find the max and indicies
  for(i=c=0; i<l; i++) {
    if(d[i] > *m) {
      // New max; store it, and reset t and c
      *m = d[i];
      t[0] = i;
      c=1;
    } else if (d[i] == *m) {
      // Tie; record the tie
      t[c++] = i;
    }
  }

  // Break any ties with a random choice
  *w = t[rnd(&R, c)];
}


// Sets m to be the max value in the size l array d, and sets w to be the index of m in d
static void Maxint(int *d, u64b l, int *m, u64b *w)
{
  int i,c,y[1024];

  // Find the max and indicies
  for(i=c=0; i<l; i++) {
    if(d[i] > *m) {
      //New max; store it, and reset t and c
      *m = d[i];
      y[0] = i;
      c=1;
    } else if (d[i] == *m) {
      // Tie; record the tie
      y[c++] = i;
    }
  }

  // Break any ties with a random choice
  *w = y[rnd(&R, c)];
}


/*
  Returns the value of the locus l in allele array c in individual i
  Individuals have 4 allele arrays:
  ->x0  // Characteristic alleles 1
  ->x1  // Characteristic alleles 2
  ->y0  // Preference alleles 1
  ->y1  // Preference alleles 2

  w stands for 'which'.  This identifies either:
  ecological traits: w == 0
  mating     traits: w == 1
  marker     traits: w == 2
  fpref      traits: w == 3
  neutral    traits: w == 4
*/
static u64b GetLocus(u64b *c, u64b l, int w)
{
  /*
  // original version below
  switch(w) {
  case 0:  return (c[l/LPW_E] & (MASK_E<<(_Ae*(l%LPW_E)))) >> (_Ae*(l%LPW_E));
  case 1:  return (c[l/LPW_P] & (MASK_P<<(_Ap*(l%LPW_P)))) >> (_Ap*(l%LPW_P));
  case 2:  return (c[l/LPW_M] & (MASK_M<<(_Am*(l%LPW_M)))) >> (_Am*(l%LPW_M));
  case 3:  return (c[l/LPW_K] & (MASK_K<<(_Ak*(l%LPW_K)))) >> (_Ak*(l%LPW_K));
  case 4:  return (c[l/LPW_F] & (MASK_F<<(_Af*(l%LPW_F)))) >> (_Af*(l%LPW_F));
  case 5:  return (c[l/LPW_N] & (MASK_N<<(_An*(l%LPW_N)))) >> (_An*(l%LPW_N));
  default: return 0;
  }*/
  switch(w) {
  case 0:  return (c[l/LPW_E] & (MASK_E<<(_Ae*(l%LPW_E)))) >> (_Ae*(l%LPW_E));
  case 1:  return (c[l/LPW_M] & (MASK_M<<(_Am*(l%LPW_M)))) >> (_Am*(l%LPW_M));
  case 2:  return (c[l/LPW_K] & (MASK_K<<(_Ak*(l%LPW_K)))) >> (_Ak*(l%LPW_K));
  case 3:  return (c[l/LPW_F] & (MASK_F<<(_Af*(l%LPW_F)))) >> (_Af*(l%LPW_F));
  case 4:  return (c[l/LPW_N] & (MASK_N<<(_An*(l%LPW_N)))) >> (_An*(l%LPW_N));
  default: return 0;
  }
}


// Returns a string n that is the same as the input string f, except that:
// %s   Is replaced by the current seed (Seed)
// %g   Is replaced by the requested number of generations (command line arg: NGen)
// %k   Is replaced by the number of characters (_k)
// %K   Is replaced by the maximum carying capacity (_K0)
// I think this is only used to generate a file name based on some 
// parameter/command-line variables.
char *FormatName(char *f)
{

  static char n[2048],b[2048];

  for(n[0]=0; *f; f++) {
    if( *f == '%' ) {
      // Special format info follows
      switch(*(++f)) {
      case 's':
        sprintf(b,"%s%llu",n,Seed);
        break;
      case 'g':
        sprintf(b,"%s%llu",n,NGen);
        break;
      case 'k':
        sprintf(b,"%s%d",n,_k);
        break;
      case 'K':
        sprintf(b,"%s%.4lf",n,((double)_K0));
        break;
      default:
        Error("Invalid name/directory format!\n");
      }
    } else {
      // No format info: just copy
      sprintf(b,"%s%c",n,*f);
    }
    // Move b back into n
    sprintf(n,"%s",b);
  }

  // n is static, so it is "safe" to return here
  return n;
}


/***************************************************************************************
 * Functions dealing with the creation and syncronization of threads
 ***************************************************************************************/

// Creates worker threads and deps, and returns when they have initialized themselves
static void CreateThreads()
{

  int i,ids[_CPUS];  pthread_t t;  pthread_attr_t a;

  // Create some condition variables and muticies for thread management
  pthread_cond_init((pthread_cond_t*)&Master,NULL);
  pthread_cond_init((pthread_cond_t*)&Slaves,NULL);
  pthread_mutex_init((pthread_mutex_t*)&Mutex,NULL);

  // Set up the thread attributes
  pthread_attr_init(&a);
  pthread_attr_setdetachstate(&a,PTHREAD_CREATE_DETACHED);
  pthread_attr_setscope(&a,PTHREAD_SCOPE_SYSTEM);

  // Create worker threads
  for(i=0; i<_CPUS; i++) {
    ids[i] = i;
    // !! It would be nice to add the thread affinity code from sched.h here, but 
    // that will only work on 2.5.8 and newer kernels...
    // WorkThreads[i] here takes it back to THREAD in main in threads.c
    if( pthread_create(&t, &a, WorkThreads[i], (void*)(ids+i)) )
      Error("Could not create thread #%d.",i);
  }

  // Wait for them to initialize
  pthread_mutex_lock((pthread_mutex_t*)&Mutex);
  for(i=0; i<_CPUS; i++) {
    while( !Done[i] )
      pthread_cond_wait((pthread_cond_t*)&Master,(pthread_mutex_t*)&Mutex);
    Done[i] = 0;
  }
  pthread_mutex_unlock((pthread_mutex_t*)&Mutex);
}


// Signals threads with flag f, and waits for them to complete
static void SignalWork(int f)
{
  
  int i;
#if _BENCHMARK
  struct timeval st,et;
  gettimeofday(&st,NULL);
  nProcessed += nCurrentIndividuals();
#endif

  // Signal worker threads and wait for them to finish
  pthread_mutex_lock((pthread_mutex_t*)&Mutex);
  for(i=0; i<_CPUS; i++)
    Go[i] = f;
  pthread_cond_broadcast((pthread_cond_t*)&Slaves);
  for(i=0; i<_CPUS; i++) {
    while( !Done[i] )
      pthread_cond_wait((pthread_cond_t*)&Master,(pthread_mutex_t*)&Mutex);
    Done[i] = 0;
  }
  pthread_mutex_unlock((pthread_mutex_t*)&Mutex);
#if _BENCHMARK
  gettimeofday(&et,NULL);
  ProcessTime += ((et.tv_sec*((u64b)1000000))+et.tv_usec) - ((st.tv_sec*((u64b)1000000))+st.tv_usec);
#endif
}


/***************************************************************************************
 * Functions which compute probibilities, properties, or prefrences of individuals
 ***************************************************************************************/

// Does any application initialization needed before starting the simulation
static void Init()
{

  u64b c,i,j,k,l;  char fn[32],n[2048];  struct timeval tv;  pIndividual t;  FILE *Parf;  void *bk;  double b;
  MaxPopSize = 0;
  
  // Do a few sanity checks
  if( _INIT_INDIVIDUALS > _MAX_INDIVIDUALS )
    Error("_INIT_INDIVIDUALS exceeds _MAX_INDIVIDUALS! Try increaseing _MAX_INDIVIDUALS.\n");
  if( ((double)_WIDTH)/_CPUS  != _WIDTH/_CPUS )
    Error("_WIDTH is not divisible by _CPUS!\n");

  // Use this to track memory useage
  bk = sbrk(0);

  // modified so that there are only 100, 010, 001, no 101, etc...
  // Fill in the Niches[] array for use by SyncIndividual() (for Log())
  for(k=0; k<_k; k++){
    //Niches[k].n = k; the line below fixes the problem with k=0 not corresponding to a niche (e.g. pm=0)
    Niches[k].n = k;
    for(i=0; i<_k; i++)
      if( (((u64b)1)<<k)&(((u64b)1)<<i) )     Niches[k].l[i] = 1;
    else                                      Niches[k].l[i] = 0;
  }

  // Set the global seed, and initialize the main thread's random number generator
  if (_RAND_SEED)
    Seed = _RAND_SEED;
  else {
    gettimeofday(&tv,NULL);
    Seed = (u64b)(tv.tv_usec|7);
  }
  initrand(&R, Seed);
  // Just go ahead and print out the seed
  printf("   Seed:           %llu\n",Seed);
  
  // Fill in the Start[][] array
  ALLOC(Start,_WIDTH*sizeof(Patch*), "Could not allocate starting patch space!\n");
  for(i=0; i<_WIDTH; i++) {
    ZALLOC(Start[i],_HEIGHT*sizeof(Patch), "Could not allocate starting patch space!\n");
    for(j=0; j<_HEIGHT; j++) {
      Start[i][j].lat = i;
      Start[i][j].lon = j;

      Start[i][j].l[0] = (byte)(0);
      Start[i][j].l[1] = (byte)(1);
      Start[i][j].l[2] = (byte)(0);
      Start[i][j].l[3] = (byte)(0);
      SyncPatch(&Start[i][j]);

      if( (i<_WIDTH) && (j<_HEIGHT) )  Start[i][j].ni = _INIT_INDIVIDUALS;

      if(Start[i][j].ni) {
#if _DISTANCEDENSITYDISPERSAL
        Start[i][j].density = ((double)(_INIT_INDIVIDUALS)) / ((double)(_K0));  // density of this patch
#endif
        ZALLOC(t, _INIT_INDIVIDUALS*sizeof(Individual), "Could not allocate starting individuals!\n");
        for(k=0; k<Start[i][j].ni; k++) {
          LinkIndividual(Start[i][j].i[k] = t+k);
          // Give this individual it's sex and other traits
          for(l=0; l<T_LENGTH; l++) Start[i][j].i[k]->d[l] = 0;
          
          
#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C
          // only in _HOMERKO, otherwise it would have initialized it to '0' as per above
          // location for fertility in _HOMERKO
          // haplosufficient
          // m locus alleles are 0: nonfunctional allele, 1: wildtype (_HOMERKO)
          for(c=0; c<_Lm; c++)    SetLocus((t+k)->m0, c, 1, 1);
          for(c=0; c<_Lm; c++)    SetLocus((t+k)->m1, c, 1, 1);
#endif
          // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING, _HOMER and _HOMERKO)
          // it is prolactin wildtype for both approaches but used differently in the two approaches
          // HOMING:                        haplosufficient fertility gene, drive is embedded in this locus
          // HOMER & HOMERKO:   haploINsufficient viability gene, drive is embedded in this locus with rescue
          for(c=0; c<_Lf; c++)    SetLocus((t+k)->f0, c, 1, 3);
          for(c=0; c<_Lf; c++)    SetLocus((t+k)->f1, c, 1, 3);

#if _SEXCHROMOSOME
// here _Ln is used as a sex chromosome; heteroz (0,1) are males; homoz. (1,1) are females.
// this is to make sure that half of pop is male, and other half is female
          if( k < (_INIT_INDIVIDUALS/2) ){
            for(c=0; c<_Ln; c++)    SetLocus((t+k)->z0, c, 0, 4);
            for(c=0; c<_Ln; c++)    SetLocus((t+k)->z1, c, 1, 4);
            (t+k)->s = 0;
          //printf("less: k: %d, x:%d \n",k, (t+k)->s);
          } else {
            for(c=0; c<_Ln; c++)    SetLocus((t+k)->z0, c, 1, 4);
            for(c=0; c<_Ln; c++)    SetLocus((t+k)->z1, c, 1, 4);
            (t+k)->s = 1;
            //printf("more: k: %d, x:%d \n",k, (t+k)->s);
          }
#else
          // if _Ln is NOT a sex chromosome, then assign sexes randomly
          (t+k)->s = rnd(&R, 2);
#endif
          (t+k)->age = 0;
        } // closes k loop
      } else {
#if _DISTANCEDENSITYDISPERSAL
        Start[i][j].density = 0.;   // density of this patch
#endif
      } // closes if(Start[i][j].ni)
    }
  }

  // Create worker threads and dependencies and wait for them to initialize
  CreateThreads();

  // Start[][] is not needed anymore, go ahead and free it
  for(i=0; i<_WIDTH; i++) {
    for(j=0; j<_HEIGHT; j++) {
      //if(Start[i][j].ei)
      //  free(Start[i][j].ei);
      free(*Start[i][j].i);
    }
    free(Start[i]);
  }
  free(Start);

  // Parse our file/dir name format strings from the parameters file
  n[0] = 0;
  sprintf(n,"%s",FormatName(_FILE_NAME));

  // Copy parameters file and append Seed there
  sprintf(fn,"cp parameters.h %s.par",n);
  system(fn);
  sprintf(fn,"%s.par",n);
  if( !(Parf=fopen(fn,"a")) ) 
    Error("Could not open file \"%s\"!\n",fn);
  fprintf(Parf,"\n//Actual Seed:           %llu\n",Seed);
  fclose(Parf);

  // Compute a value for _b which will allow view to scale properly
  //b = 1.0;  //original line
  b = _b/2.;

#if _SYS
  // Open sys file
  sprintf(fn,"%s.sys",n);
  if( !(Sysf=fopen(fn,"w")) ) 
    Error("Could not open file \"%s\"!\n",fn);
  // Write the header
  fprintf(Sysf,"%d %lf %lf %llu %llu\n",_k,((double)_K0),b,(u64b)(_WIDTH),(u64b)(_HEIGHT));
#endif
#if _DEME
  // Open deme file
  sprintf(fn,"%s.deme",n);
  if( !(Demef=fopen(fn,"w")) ) 
    Error("Could not open file \"%s\"!\n",fn);
  // Write the header
  fprintf(Demef,"%d %lf %lf %llu %llu\n",_k,((double)_K0),b,(u64b)(_WIDTH),(u64b)(_HEIGHT));
#endif
#if _DISTANCE
  // Open distance file
  sprintf(fn,"%s.distance",n);
  if( !(Distancef=fopen(fn,"w")) )
    Error("Could not open file \"%s\"!\n",fn);
  fprintf(Distancef,"%d\n", ((int)((((double)(_DRANGE))-1.)/2.)) + 1);
  fprintf(Distancef, "gen: %d ", 0);
  // this part is for viewing
  // otherwise since dispersal occurs at gen 0, it is viewed at gen 0
  // we want to see it 'after' it happened at gen 1
  int ii, jj;
  for( ii=0; ii<10; ii++){
    for( jj=0; jj<((int)((((double)(_DRANGE))-1.)/2.)) + 1; jj++){
      fprintf(Distancef, "0 ");
    }
  }
  fprintf(Distancef, "\n");
#endif
#if _IND
  // Open individual file
  sprintf(fn,"%s.ind",n);
  if( !(Indf=fopen(fn,"w")) ) 
    Error("Could not open file \"%s\"!\n",fn);
  // Write the header
  // fprintf(Indf,"%d %lf %lf %llu %llu\n",_k,((double)_K0),b,(u64b)(_WIDTH),(u64b)(_HEIGHT));
#endif

  // How much memory are we using now?  
  printf("   Mem:            %lluK\n",(((u64b)sbrk(0))-((u64b)bk))/(1<<10));
}


#if (_SYS || _DEME )
// Writes system level data to a ".sys" and a ".deme" file
static void LogSysDeme()
{
  u64b   c,i,j,k,l,max_i,max_gen,sys=1;
  u64b   genotypepop[8];
  double gcc, amc, amk;
  double cc=0.,pr=0., match[2], X, x[_k], x2[_k], f[_k],t=0.,mt=0.,kt=0.,ft=0.;
  pPatch p;
  int    lat, lon;

  // Every generation in the log files is now preceded with the tag 'gen: '.
  // This should help to solve some of the problems with data file "skewing".
#if _SYS
  if( !(Gen%_SYS ) )
    fprintf(Sysf, "gen: ");
#endif
#if _DEME
  if( !(Gen%_DEME) )
    fprintf(Demef, "gen: ");
#endif

#if _SYS
  if( !(Gen%_SYS ) ){
    /* If there are no surviving individuals, just print a header and move on */
    // ReadSys() in view.c doesn't like moving on...
    if(!nCurrentIndividuals()) {
      fprintf(Sysf, "%llu 0 \n", Gen);
      //Sys_logs++; // doesn't like incrementing when there is extinction. so leave it silenced here
      sys = 0;
    } else {
      // Init generation level vars
      amc = 0.0;
      memset(genotypepop, 0, sizeof(genotypepop));
    }
  }
#endif

  //int ttt;
  int pall, ptotal;
  ptotal = 8;
  // mean for all x traits
  
  
  for(lat=0; lat<WIDTH; lat++){
    for(lon=0; lon<HEIGHT; lon++){
      p = &Threads[0]->Current[lat][lon];
      //printf("lat: %d lon: %d\n", lat, lon);
      for(pall=0; pall<ptotal; pall++){
        p->AllSpecies[pall] = 0;
        p->speciespops[pall] = 0;
      }
      
      int pop = 0;
      // loop through individuals
      for(k=0; k<p->ni; k++){
        //printf("k:%d (%d, %d) pop: %d\n", k, lat, lon, pop);
#if _HOMING_A
        // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING)
        // if it has allele '2' it has drive
        if( *(p->i[k]->f0) == 2 || *(p->i[k]->f1) == 2 ){
            p->AllSpecies[1] = 1;
            p->speciespops[1]++;
            pop++;
        } else if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3 ){
            p->AllSpecies[3] = 1;
            p->speciespops[3]++;
            pop++;
        } else if( *(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4 ){
            p->AllSpecies[4] = 1;
            p->speciespops[4]++;
            pop++;
        }
#endif
        //  printf("in between\n");
#if _HOMING_B
        // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING)
        // if male
        //printf("*k:%d (%d, %d) pop: %d sex: %d\n", k, lat, lon, pop, p->i[k]->s);
        //printf("%d\n", p->i[k]->s);
        if( p->i[k]->s == 0 ){
          //printf("here1?\n");
          // if it has allele '2' it has drive
          if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[5] = 1;
            p->speciespops[5]++;
            pop++;
          } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
            p->AllSpecies[4] = 1;
            p->speciespops[4]++;
            pop++;
          } else if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ||
                     (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4) ||
                     (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 4) ||
                     (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[6] = 1;
            p->speciespops[6]++;
            pop++;
          } else {
            p->AllSpecies[7] = 1;
            p->speciespops[7]++;
            pop++;
          }
        } else {
          //printf("here2\n");
          // if it has allele '2' it has drive
          if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[1] = 1;
            p->speciespops[1]++;
            pop++;
          } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
            p->AllSpecies[0] = 1;
            p->speciespops[0]++;
            pop++;
          } else {
            p->AllSpecies[3] = 1;
            p->speciespops[3]++;
            pop++;
          }
        }
        //printf("here end?\n");
#endif
          
#if _HOMING_C
        // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING)
        // if male
        //printf("*k:%d (%d, %d) pop: %d sex: %d\n", k, lat, lon, pop, p->i[k]->s);
        //printf("%d\n", p->i[k]->s);
        if( p->i[k]->s == 0 ){
          //printf("here1?\n");
          // if it has allele '2' it has drive
          if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[5] = 1;
            p->speciespops[5]++;
            pop++;
          } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
            p->AllSpecies[4] = 1;
            p->speciespops[4]++;
            pop++;
          } else if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ||
                     (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4) ||
                     (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 4) ||
                     (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[6] = 1;
            p->speciespops[6]++;
            pop++;
          } else {
            p->AllSpecies[7] = 1;
            p->speciespops[7]++;
            pop++;
          }
        } else {
          //printf("here2\n");
          // if it has allele '2' it has drive
          if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[1] = 1;
            p->speciespops[1]++;
            pop++;
          } else if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ||
                     (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4) ||
                     (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 4) ||
                     (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[2] = 1;
            p->speciespops[2]++;
            pop++;
          } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
            p->AllSpecies[0] = 1;
            p->speciespops[0]++;
            pop++;
          } else {
            p->AllSpecies[3] = 1;
            p->speciespops[3]++;
            pop++;
          }
        }
        //printf("here end?\n");
#endif
          
#if _HOMER
        // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMER)
        //printf("*k:%d (%d, %d) pop: %d sex: %d\n", k, lat, lon, pop, p->i[k]->s);
        //printf("%d\n", p->i[k]->s);
        // if it has allele '2' it has drive
        if( *(p->i[k]->f0) == 2 || *(p->i[k]->f1) == 2 ){
            p->AllSpecies[1] = 1;
            p->speciespops[1]++;
            pop++;
        } else if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3 ){
            p->AllSpecies[3] = 1;
            p->speciespops[3]++;
            pop++;
        } else if( *(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4 ){
            p->AllSpecies[4] = 1;
            p->speciespops[4]++;
            pop++;
        }
        //printf("here end?\n");
#endif

#if _HOMER_HS
        // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMER)
        // if male
        //printf("*k:%d (%d, %d) pop: %d sex: %d\n", k, lat, lon, pop, p->i[k]->s);
        //printf("%d\n", p->i[k]->s);
        if( p->i[k]->s == 0 ){
          //printf("here1?\n");
          // if it has allele '2' it has drive
          if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[5] = 1;
            p->speciespops[5]++;
            pop++;
          } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
            p->AllSpecies[4] = 1;
            p->speciespops[4]++;
            pop++;
          } else if( *(p->i[k]->f0) == 5 || *(p->i[k]->f1) == 5 ){
            p->AllSpecies[7] = 1;
            p->speciespops[7]++;
            pop++;
          } else {
            p->AllSpecies[6] = 1;
            p->speciespops[6]++;
            pop++;
          }
        } else {
          //printf("here2\n");
          // if it has allele '2' it has drive
          if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
            p->AllSpecies[1] = 1;
            p->speciespops[1]++;
            pop++;
          } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
            p->AllSpecies[0] = 1;
            p->speciespops[0]++;
            pop++;
          } else if( *(p->i[k]->f0) == 5 || *(p->i[k]->f1) == 5 ){
            p->AllSpecies[3] = 1;
            p->speciespops[3]++;
            pop++;
          } else {
            p->AllSpecies[2] = 1;
            p->speciespops[2]++;
            pop++;
          }
        }
        //printf("here end?\n");
#endif

          
#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C
        // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
        // m locus alleles are 0: nonfunctional, 1: wildtype, 1: drive
        // anything with '4' in f is inviable so don't count
        // '11' in f is default so don't count (for deme colors)
        /*
        if( *(p->i[k]->f0) == 2 || *(p->i[k]->f1) == 2 ){         // at least one copy of drive
          if( *(p->i[k]->m0) == 1 || *(p->i[k]->m1) == 1  ){      // at least one functional fertility
            p->AllSpecies[1] = 1;
            p->speciespops[1]++;
            pop++;
            continue;
          } else {
            // infertile
            p->AllSpecies[2] = 1;
            p->speciespops[2]++;
            pop++;
            continue;
          }
        } else if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3 ){  // at least one copy of cut functional (resistant)
          if( *(p->i[k]->m0) == 1 || *(p->i[k]->m1) == 1  ){      // at least one functional fertility  p->AllSpecies[3] = 1;
            p->AllSpecies[3] = 1;
            p->speciespops[3]++;
            pop++;
            continue;
          } else {
            // infertile
            p->AllSpecies[4] = 1;
            p->speciespops[4]++;
            pop++;
            continue;
          }
        } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){  // at least one copy of cut functional (resistant)
          if( *(p->i[k]->m0) == 0 || *(p->i[k]->m1) == 0  ){      // at least one functional fertility  p->AllSpecies[3] = 1;
            // infertile
            p->AllSpecies[5] = 1;
            p->speciespops[5]++;
            pop++;
            continue;
          }
        }*/
        if( *(p->i[k]->f0) == 2 || *(p->i[k]->f1) == 2 ){         // at least one copy of drive
          p->AllSpecies[1] = 1;
          p->speciespops[1]++;
          pop++;
          continue;
        } else if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3 ){  // at least one copy of cut functional (resistant)
          p->AllSpecies[3] = 1;
          p->speciespops[3]++;
          pop++;
          continue;
        } else if( *(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4 ){  // at least one copy of cut functional (resistant)
          p->AllSpecies[5] = 1;
          p->speciespops[5]++;
          pop++;
          continue;
        }
#endif
        
      }
      //printf("here3?\n");
      //if( pop > p->ni ) printf("\nno. of inds assigned to a species is: %d, but N in patch is: %llu", pop, p->ni);

      if( pop == 0 ){
        p->SpeciesID = 0;
      } else {
        // set the species id of the patch to the most abundant species
        int m = 0;
        u64b wh;
        Maxint(p->speciespops, ptotal, &m, &wh);
        p->SpeciesID = wh;
        //if( wh != 0 && wh != 4 ) printf("(%d, %d) species id is set to: %d (%d/%d)\n\n", lat, lon, wh, p->speciespops[wh], p->ni);
      }
      //printf("here4?\n");
    }
    //printf("here5?\n");
  }
  //printf("here6?\n");

  for(gcc=0.,amc = 0.,amk=0.,max_i=max_gen=c=0; c<_CPUS; c++){
    for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
        p = &Threads[c]->Current[i][j];

#if _DEME
        // Init deme-level vars
        if( !(Gen%_DEME) ) {
          cc = 0.;
          pr = 0.;
          memset(x,     0, sizeof(x));
          memset(x2,    0, sizeof(x2));
          mt = kt = ft = 0.0;
          memset(match, 0, sizeof(match));
        }
#endif
#if _SYS
        if( !(Gen%_SYS) && sys ) {
          // Find max individual count per patch
          if( p->ni > max_i )
            max_i = p->ni;
        }
#endif
        for(k=0; k<p->ni; k++) {
#if _DEME
          if( !(Gen%_DEME) ){
            //cc += K(p->i[k],p);
            //pr += P(p->i[k],p);
            // Sum over all matching fit/mig
            // NOT relevant
            //match[0] += Qf(p->i[k],p);
            //match[1] += Qm(p->i[k],p);
            match[0] += 1;
            match[1] += 1;

            // Sum over all traits
            for(l=0; l<_k; l++) {
              X = p->i[k]->x[l];
              //Y = p->i[k]->y[l];
              x[l] += X;
              //y[l] += Y;
              x2[l] += X*X;
              //y2[l] += Y*Y;
            }

#if _HOMING_A || _HOMING_B || _HOMING_C || _HOMER || _HOMER_HS
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
            // drive frequency
            if( *(p->i[k]->f0) == 2){
              gcc++;
            }
            if( *(p->i[k]->f1) == 2){
              gcc++;
            }
            // wildtype frequency
            if( *(p->i[k]->f0) == 1){
              amc++;
            }
            if( *(p->i[k]->f1) == 1){
              amc++;
            }
#if _EXPOSURE
            // no-rescue frequency
            if( *(p->i[k]->f0) == 5 ){
              amk++;
            }
            if( *(p->i[k]->f1) == 5 ){
              amk++;
            }
#else
            // resistance frequency (but 3 is functional if _RESISTANCE)
            if( *(p->i[k]->f0) == 3){
              amk++;
            }
            if( *(p->i[k]->f1) == 3){
              amk++;
            }
            // nonfunctional resistant essential gene
            if( *(p->i[k]->f0) == 4){
              amk++;
            }
            if( *(p->i[k]->f1) == 4){
              amk++;
            }
#endif
#endif
            
#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C
            // m locus alleles are 0: nonfunctional, 1: wildtype (_HOMERKO)
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
            // infertile allele frequency in 'm'
            if( *(p->i[k]->f0) == 1){
              amc++;
            }
            if( *(p->i[k]->f1) == 1){
              amc++;
            }
            // drive frequency in 'f'
            if( *(p->i[k]->f0) == 2){
              gcc++;
            }
            if( *(p->i[k]->f1) == 2){
              gcc++;
            }
            // infertile
            if( *(p->i[k]->m0) == 0){
              amk++;
            }
            if( *(p->i[k]->m1) == 0){
              amk++;
            }
#endif
            
#if _HOMING_A
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING)
            // if female
            if( p->i[k]->s == 1 ){
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
                genotypepop[1]++;
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[0]++;
              } else if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ||
                         (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4) ||
                         (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 4) ||
                         (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 2) ){
                genotypepop[2]++;
              } else {
                genotypepop[3]++;
              }
            } else {
              // if male
              genotypepop[6] = 0; // no infertile males in 1A
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
                genotypepop[5]++;
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[4]++;
              } else {
                genotypepop[7]++;
              }
            }
#endif

#if _HOMING_B
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING)
            // if male
            if( p->i[k]->s == 0 ){
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
                genotypepop[5]++;
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[4]++;
              } else if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ||
                         (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4) ||
                         (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 4) ||
                         (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 2) ){
                genotypepop[6]++;
              } else {
                genotypepop[7]++;
              }
            } else {
              // if female
              genotypepop[2] = 0; // no infertile females in 1A
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
                genotypepop[1]++;
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[0]++;
              } else {
                genotypepop[3]++;
              }
            }
#endif

#if _HOMING_C
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING)
            // if male
            if( p->i[k]->s == 0 ){
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
                genotypepop[5]++;
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[4]++;
              } else if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ||
                         (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4) ||
                         (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 4) ||
                         (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 2) ){
                genotypepop[6]++;
              } else {
                genotypepop[7]++;
              }
            } else {
              // if female
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ||
                    (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 4) ||
                    (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 4) ||
                    (*(p->i[k]->f0) == 4 && *(p->i[k]->f1) == 2) ){
                  genotypepop[2]++; // no infertile females in 1A
                  // if it has allele '2' it has drive
              } else if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) ){
                genotypepop[1]++;
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[0]++;
              } else {
                genotypepop[3]++;
              }
            }
#endif


#if _HOMER
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING and _HOMER)
            // if male
            if( p->i[k]->s == 0 ){
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
                genotypepop[5]++; // drive
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[4]++; // wildtype
              } else {
                genotypepop[6]++; // others
                // if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3) printf("m:%llu/%llu\n", *(p->i[k]->f0), *(p->i[k]->f1));
              }
            } else {
              // if female
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
                genotypepop[1]++; // drive
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[0]++; // wildtype
              } else {
                genotypepop[2]++; // others
                // if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3) printf("f:%llu/%llu\n", *(p->i[k]->f0), *(p->i[k]->f1));
              }
            }
#endif

              
#if _HOMER_HS
#if _EXPOSURE
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING and _HOMER)
            // if male
            if( p->i[k]->s == 0 ){
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
                genotypepop[5]++; // drive
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[4]++; // wildtype
              } else if( *(p->i[k]->f0) == 5 || *(p->i[k]->f1) == 5 ){
                genotypepop[7]++; // no rescue
              } else {
                genotypepop[6]++; // others
                // if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3) printf("m:%llu/%llu\n", *(p->i[k]->f0), *(p->i[k]->f1));
              }
            } else {
              // if female
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
                genotypepop[1]++; // drive
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[0]++; // wildtype
              } else if( *(p->i[k]->f0) == 5 || *(p->i[k]->f1) == 5 ){
                genotypepop[3]++; // no rescue
              } else {
                genotypepop[2]++; // others
                // if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3) printf("f:%llu/%llu\n", *(p->i[k]->f0), *(p->i[k]->f1));
              }
            }
#else
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant (_HOMING and _HOMER)
            // if male
            if( p->i[k]->s == 0 ){
              // if it has allele '2' it has drive
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
                genotypepop[5]++; // drive
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[4]++; // wildtype
              } else {
                genotypepop[6]++; // others
                // if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3) printf("m:%llu/%llu\n", *(p->i[k]->f0), *(p->i[k]->f1));
              }
            } else {
              // if female
              if( (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 1) || (*(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 2) || (*(p->i[k]->f0) == 2 && *(p->i[k]->f1) == 2) ){
                genotypepop[1]++; // drive
              } else if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
                genotypepop[0]++; // wildtype
              } else {
                genotypepop[2]++; // others
                // if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3) printf("f:%llu/%llu\n", *(p->i[k]->f0), *(p->i[k]->f1));
              }
            }
#endif
#endif


#if _HOMERKO_A || _HOMERKO_B || _HOMERKO_C
            // m locus alleles are 0: nonfunctional, 1: wildtype (_HOMERKO)
            // f locus alleles are 1: wildtype, 2: drive, 3: functional resistant, 4: nonfunctional resistant
            /*
            if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
              if( p->i[k]->s == 1) genotypepop[0]++;   //wildtype at least fertile female
              if( p->i[k]->s == 0) genotypepop[4]++;
            } else if( *(p->i[k]->f0) == 2 || *(p->i[k]->f1) == 2 ){ // if it has allele '2' it has drive
              if( p->i[k]->s == 1) genotypepop[1]++;   //drive and fertile female
              if( p->i[k]->s == 0) genotypepop[5]++;   //drive and fertile male
            } else if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3 ){ // if it has allele '3' it has r1
              if( p->i[k]->s == 1) genotypepop[2]++;   //r1
              if( p->i[k]->s == 0) genotypepop[6]++;
            } else if( *(p->i[k]->f0) == 4 || *(p->i[k]->f1) == 4 ){ // if it has allele '4' it has r2
              if( p->i[k]->s == 1) genotypepop[3]++;   //r2
              if( p->i[k]->s == 0) genotypepop[7]++;
            }
            */
              
            if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
              if( p->i[k]->s == 1) genotypepop[0]++;   //wildtype at least fertile female
              if( p->i[k]->s == 0) genotypepop[4]++;
            } else if( *(p->i[k]->f0) == 2 || *(p->i[k]->f1) == 2 ){ // if it has allele '2' it has drive
              if( p->i[k]->s == 1 ){ // drive and female
                if( *(p->i[k]->m0) == 1 || *(p->i[k]->m1) == 1 ){
                  genotypepop[1]++;  //drive and fertile female
                } else {
                  genotypepop[3]++;  //drive and infertile female
                }
              }
              if( p->i[k]->s == 0) genotypepop[5]++;   //drive and fertile male
            } else if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3 ){ // if it has allele '3' it has r1
              if( p->i[k]->s == 1) genotypepop[2]++;   //r1 female
              if( p->i[k]->s == 0) genotypepop[6]++;   //r1 male
            } else { // if it has allele '4' it has r2 both sexes
              genotypepop[7]++;
            }
            
            /*
            // this is to get detailed fertile female numbers
            if( *(p->i[k]->f0) == 1 && *(p->i[k]->f1) == 1 ){
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 1 || *(p->i[k]->m1) == 1 ) ) genotypepop[0]++;   //wildtype at least fertile female
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 0 && *(p->i[k]->m1) == 0 ) ) genotypepop[4]++;
            } else if( *(p->i[k]->f0) == 2 || *(p->i[k]->f1) == 2 ){ // if it has allele '2' it has drive
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 1 || *(p->i[k]->m1) == 1 ) ) genotypepop[1]++;   //drive and fertile female
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 0 && *(p->i[k]->m1) == 0 ) ) genotypepop[5]++;   //drive and fertile male
            } else if( *(p->i[k]->f0) == 3 || *(p->i[k]->f1) == 3 ){ // if it has allele '3' it has r1
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 1 || *(p->i[k]->m1) == 1 ) ) genotypepop[2]++;   //r1
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 0 && *(p->i[k]->m1) == 0 ) ) genotypepop[6]++;
            } else if( *(p->i[k]->f0) == 4 || *(p->i[k]->f1) == 4 ){ // if it has allele '4' it has r2
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 1 || *(p->i[k]->m1) == 1 ) ) genotypepop[3]++;   //r2
              if( p->i[k]->s == 1 && ( *(p->i[k]->m0) == 0 && *(p->i[k]->m1) == 0 ) ) genotypepop[7]++;
            } */
#endif
            
            mt  += p->i[k]->m;
            kt  += p->i[k]->k;
            ft  += p->i[k]->f;
          }
#endif

        }

#if _SYS
        if( !(Gen%_SYS) && sys ) {
          // Find max genotype count to scale
          for(l=0; l<8; l++) {
            if( genotypepop[l] > max_gen ) max_gen = genotypepop[l];
          }
        }
#endif

#if _DEME
        if( !(Gen%_DEME) ) {
          // Find most prefered niche
          if ( (k = p->ni) ){ t = 0; Max(f, NICHES, &t, &l); t = l; }
          // Divide by nIndividuals for average and find variances
          k        = p->ni;

          cc      /= k;
          pr      /= k;
          match[0]/= k;
          match[1]/= k;

          for(l=0; l<_k; l++) {
            x[l]  /= k;
            x2[l]  = x2[l]/k - x[l]*x[l];
          }
          mt /= k;
          kt /= k;
          ft /= k;
#endif

#if _DEME
          // Write this patch's stats
          fprintf(Demef, "%llu %llu %llu %lf %lf %lf %lf %llu ", Gen, p->ni, ((u64b)t), cc, pr, match[0], match[1], p->n);
          // Write niche type
          for(l=0; l<_k; l++) fprintf(Demef, "%d ", p->l[l]);
          // Write average traits and variances
          for(l=0; l<_k; l++) fprintf(Demef, "%lf ", x[l]);
          for(l=0; l<_k; l++) fprintf(Demef, "%lf ", x2[l]);
          // In an attempt to make the viewer compatible with both sexual versions
          // as well as the non-sexual case, I will embed a flag for each option.
          fprintf(Demef, "1 %lf %lf ", mt, kt);
          fprintf(Demef, "1 %lf ", ft);
          fprintf(Demef, "1 ");
          // other info could be present, but if this is 0, NOT a species!!!
          fprintf(Demef, "%llu ", p->SpeciesID);
          // Terminate with a newline
          fprintf(Demef, "\n");
        }
#endif
      }
    }
  }

#if _SYS
  if( !(Gen%_SYS) && sys ) {

    // Divide by nIndividuals() to find averages
    i = nCurrentIndividuals();
    
    gcc            /= ((double)(2*i));
    amc            /= ((double)(2*i));
    amk            /= ((double)(2*i));

#if _HOMER
    if( gen_n == 0 && gcc >= 0.90 ){
      gen_n = Gen;
    }
    if( gen_nf == 0 && gcc >= 0.95 ){
      gen_nf = Gen;
    }
    if( gen_nn == 0 && gcc >= 0.99 ){
      gen_nn = Gen;
    }
#endif
#if _EXPOSURE
      if( gen_exp == 0 && gcc >= _ExposureFrequency ){
        gen_exp = Gen;
      }
#endif
      
    fprintf(Sysf, "%llu %llu %lf %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu %lf %lf %d %d %d %d",
            Gen, i, gcc, max_gen, max_i, genotypepop[0], genotypepop[1], genotypepop[2], genotypepop[3], genotypepop[4], genotypepop[5], genotypepop[6], genotypepop[7], amc, amk, OffspringPopSize, Mortality, MultiplePaternity, SinglePaternity);
    // Terminate with a newline
    fprintf(Sysf,"\n");
    fflush(Sysf);
    // Incriment counter
    Sys_logs++;
  }
#endif
#if _DEME
  // Increment counter
  if( !(Gen%_DEME) ) {
    fflush(Demef);
    Deme_logs++;
  }
#endif


}
#endif


#if _IND
// Dumps individual data to a ".ind" file
static void LogInd()
{
  u64b t,c,i,j,k,n,pf,pm,x,y;  pPatch p;  pIndividual u;  double f,m;
  u64b l;

  if(!nCurrentIndividuals()) {
    fprintf(Indf, "%llu 0\n", Gen);
    Ind_logs++;
  } else {
  //if( !(_IND_SAMPLE) || (_IND_SAMPLE < nCurrentIndividuals()) ) {
    //////////////////////////////////////////////////////////////////////
    // Log normally without any subsample
    //////////////////////////////////////////////////////////////////////
    //for(c=0; c<_CPUS; c++)
      for(i=0; i<WIDTH; i++) {
        for(j=0; j<HEIGHT; j++) {
          p = &Threads[0]->Current[i][j];
          for(k=0; k<p->ni; k++) {
            // Print allelic values
            fprintf(Indf, "%llu %llu %llu 1 ", Gen, i, j);
            // Print trait values
            for(l=0; l<_k; l++) fprintf(Indf,"%lf ",p->i[k]->x[l]);
            // print the loci
            for(l=0; l<_k; l++) for(n=0; n<_Le; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->x0,(l*_Le)+n,0));
            for(l=0; l<_k; l++) for(n=0; n<_Le; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->x1,(l*_Le)+n,0));
            // Print trait values
            fprintf(Indf,"%lf ",p->i[k]->m);
            // print the loci
            for(n=0; n<_Lm; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->m0,n,2));
            for(n=0; n<_Lm; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->m1,n,2));
            // Print trait values
            fprintf(Indf,"%lf ",p->i[k]->k);
            // print the loci
            for(n=0; n<_Lk; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->k0,n,3));
            for(n=0; n<_Lk; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->k1,n,3));
            // Print trait values
            fprintf(Indf,"%lf ",p->i[k]->f);
            // print the loci
            for(n=0; n<_Lf; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->f0,n,4));
            for(n=0; n<_Lf; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->f1,n,4));
            // Neutral traits only have one character
            for(n=0; n<_Ln; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->z0,n,5));
            for(n=0; n<_Ln; n++) fprintf(Indf,"%llu ",GetLocus(p->i[k]->z1,n,5));
            // Print sex
            fprintf(Indf,"%llu ",p->i[k]->s);

            fprintf(Indf,"%llu ",p->i[k]->age);

            fprintf(Indf, "\n");
            fflush(Indf);
          }
        }
      }
  // Incriment counter
  Ind_logs++;
  }
}
#endif



// Writes any necessary data to a log file(s)
static void Log()
{
  static u32b save[59];
  // int n;
  // for(n=0; n<59; n++) printf("%llu ", R.rtab[n]);
  // printf(" - before save\n");
  sv_rnd( &R, save );
    
  // Write to file(s)
#if (_SYS || _DEME)
  if( (_SYS && !(Gen%_SYS)) || (_DEME && !(Gen%_DEME)) ) LogSysDeme();
#endif
#if _IND
  if( !(Gen%_IND) )       LogInd();
#endif
  

  int year;
  year = ((int) ((double)(Gen))/((double)(_MATINGCYCLES)) );

  if( nCurrentIndividuals() > MaxPopSize ){ MaxPopSize = nCurrentIndividuals(); }
  // Messages for stdout/stderr
  printf("   cycle = %llu        year = %d        p = %llu \n", Gen, year, nCurrentIndividuals());
  rst_rnd( &R, save );
  //int n;
  //for(n=0; n<59; n++) printf("%llu ", R.rtab[n]);
  //  printf(" - after restore\n\n\n");
}

// Advances the simulation one generation
static void Generation()
{

  // Tell worker threads to start working
  SignalWork(1);

  // Advance generation counter
  Gen++;
}


/***************************************************************************************
 * Main and signal handler...
 ***************************************************************************************/

// Function to be called when the simulation exits
static void Cleanup()
{

  // Print out summary and benchmark data if needed
  printf("* Done\n");
  printf("\n* Summary:\n");
//  printf("   Niche version:  %s\n",VERSION);

  printf("   Ran for:                   %llu mating cycles\n", Gen);
  printf("   Seed:                      %llu\n", Seed);
  printf("   MaxPopSize:                %d\n", MaxPopSize);
  printf("   Mating mothers with drive: %d\n", MatingMotherwithDrive);
  printf("   Mating fathers with drive: %d\n", MatingFatherwithDrive);
      
#if _HOMER
  printf("   T(90):          %d\n",gen_n);
  printf("   T(95):          %d\n",gen_nf);
  printf("   T(99):          %d\n",gen_nn);
#endif
#if _EXPOSURE
  printf("   T(exp):         %d\n",gen_exp);
#endif
#if _BENCHMARK
  printf("\n* Benchmark:\n   %lf Seconds spent processing %llu individuals.\n   %lf Seconds per million.\n",
          ProcessTime/1000000.,nProcessed,((double)ProcessTime)/nProcessed);
#endif

  // Close the log file(s)
#if _SYS
  fprintf(Sysf,"\n%llu", Sys_logs);
  fclose(Sysf);
#endif
#if _DISTANCE
  fprintf(Distancef,"\n%llu", Distance_logs+1);
  fclose(Distancef);
#endif
#if _DEME
  fprintf(Demef,"\n%llu", Deme_logs);
  fclose(Demef);
#endif
#if _IND
  fprintf(Indf,"\n%llu", ((u64b)(NGen/_IND))+1);
  fclose(Indf);
#endif
  // All allocated memory and threads will be released by the OS upon process termination,
  // so we don't need to bother with freeing or killing threads here
}


// Signal handler to catch SIGTERM and SIGINT
static void Signals(int arg)
{

  volatile static int called=0;

  if(!called) {
    called = 1;
    fprintf(stderr, "* Signal received; exiting gracefully.\n");
    Cleanup();
    exit(0);
  }
}


// Application entry point
int main(int argc, char **argv)
{
  FILE *d;
  
  // Check command line args
  if(argc < 2) {
#if _ARBITER
    NGen = _GENERATIONS;
#else
    Error("useage:\n\t%s <ngen>\n\nwhere <ngen> is the number of generations to simulate.\n",*argv);
#endif
  } else {
    sscanf(argv[1],"%llu",&NGen);
  }

#if _ARBITER
  arbiter_init();
#endif

  // Install signal handler
  signal(SIGTERM, Signals);
  signal(SIGINT,  Signals);

  // Initialize and run the simulation
  printf("* Initializing...\n");
  // printf("   Niche version:  %s\n",VERSION);
  

  printf("   Running for:    %llu mating cycles \n",NGen);
  //printf("%lf\n",_LOFPROB);
  Init();
  
  //int n;
  //for(n=0; n<59; n++) printf("%llu ", R.rtab[n]);
  //  printf(" - after restore\n\n\n");
  
  printf("* Done.\n\n* Running simulation:\n");
  // Log() should be after LoadCheckpoint() since it shouldn't write to the files as if gen=0 if there is already an indall file
  // but it shouldn't also write the same generation twice, i.e. gen it resumes, and before advances with Generation() below
  // printf("this is before the FIRST Log()\n");
  
  Log();
  
  while( (Gen < NGen) && nCurrentIndividuals() ) {
    Generation();

#if _ARBITER
    char string[128];
    if( !(Gen%_SYS) ){
      sprintf(string,"g-%llu\n",Gen);
      arbiter_custom(string);
    }
#endif
    
    Log();
  }
  
  if( !nCurrentIndividuals() )
    printf("* Extinction.\n");

  // Clean up and return success
  Cleanup();
  d = fopen( ".done", "wt" );
  fclose(d);

#if _ARBITER
  arbiter_finished();
#endif
  return 0;
}

#endif
