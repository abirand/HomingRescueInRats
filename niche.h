#ifndef __NICHE_H
#define __NICHE_H

#include "parameters.h"
#include "types.h"

#define NICHE_LITTLE_ENDIAN 1
#define NICHE_BIG_ENDIAN    2

#define MALE                0
#define FEMALE              1

// _WIDTH and _HEIGHT as seen from an individual thread's perspective
#define WIDTH               (_WIDTH/_CPUS)
#define HEIGHT              _HEIGHT

// Maximum population sizes
// #define MAX_POP          (_WIDTH*_HEIGHT*_MAX_INDIVIDUALS)
#define MAX_THREAD_POP      (_WIDTH*_HEIGHT*_MAX_INDIVIDUALS)

// Number of different niches possible
#define NICHES              _k

// Number of u64b's needed per half characteristic-set
#define LENGTH_E            (((_k*_Le*_Ae)/64)+((((_k*_Le*_Ae)%64)>0)*1))
#define LENGTH_M            (((_Lm*_Am)/64)+((((_Lm*_Am)%64)>0)*1))
#define LENGTH_K            (((_Lk*_Ak)/64)+((((_Lk*_Ak)%64)>0)*1))
#define LENGTH_F            (((_Lf*_Af)/64)+((((_Lf*_Af)%64)>0)*1))
#define LENGTH_N            (((_Ln*_An)/64)+((((_Ln*_An)%64)>0)*1))
#define T_LENGTH            ((LENGTH_E<<1)+(LENGTH_M<<1)+(LENGTH_K<<1)+(LENGTH_F<<1)+(LENGTH_N<<1))

// Number of bytes that actually need to be processed
#define PLENGTH_E           (((_k*_Le*_Ae)/8)+((((_k*_Le*_Ae)%8)>0)*1))
#define PLENGTH_M           (((_Lm*_Am)/8)+((((_Lm*_Am)%8)>0)*1))
#define PLENGTH_K           (((_Lk*_Ak)/8)+((((_Lk*_Ak)%8)>0)*1))
#define PLENGTH_F           (((_Lf*_Af)/8)+((((_Lf*_Af)%8)>0)*1))
#define PLENGTH_N           (((_Ln*_An)/8)+((((_Ln*_An)%8)>0)*1))

// Number of bits in a u64b
#define BITS                (8*sizeof(u64b))

// Loci per byte; used by fast-recombination
#define FIELDS_E            (8/_Ae)
#define FIELDS_M            (8/_Am)
#define FIELDS_K            (8/_Ak)
#define FIELDS_F            (8/_Af)
#define FIELDS_N            (8/_An)

// Mask for pulling out a locus, and "loci per u64b(word)"
#define MASK_E              ((((u64b)1)<<_Ae)-1)
#define MASK_M              ((((u64b)1)<<_Am)-1)
#define MASK_K              ((((u64b)1)<<_Ak)-1)
#define MASK_F              ((((u64b)1)<<_Af)-1)
#define MASK_N              ((((u64b)1)<<_An)-1)
#define LPW_E               (BITS/_Ae)
#define LPW_M               (BITS/_Am)
#define LPW_K               (BITS/_Ak)
#define LPW_F               (BITS/_Af)
#define LPW_N               (BITS/_An)

// Nifty little macros to make malloc a little cleaner
#define ALLOC(t,s,e)        { if( !((t)=malloc(s)) )    Error(e);                         }
#define ZALLOC(t,s,e)       { if( !((t)=malloc(s)) )    Error(e); else memset((t),0,(s)); }
#define REALLOC(t,s,e)      { if( !((t)=realloc(t,s)) ) Error(e);                         }

typedef  void* (*WorkThread)(void*);

// Structure to represent an individual
typedef struct {
  u64b      d[T_LENGTH];                  // Space for x, y, m, k, and z below
  u64b     *x0;                           // d; ecological characters
  u64b     *x1;                           // d; ecological characters
  double    x[_k];                        // Scaled characteristis
  u64b     *m0;                           // d; gene-drive: drive gene (old: mating (choosiness in females))
  u64b     *m1;                           // d; gene-drive: drive gene (old: mating (choosiness in females))
  u64b     *k0;                           // d; marker in males (in females only if similarity based mating)
  u64b     *k1;                           // d; marker in males (in females only if similarity based mating)
  double    m;                            // Scaled drive (mating (choosiness in females))
  double    k;                            // Scaled marker (in males, in females only if no sexual selection)
  u64b     *f0;                           // d; gene-drive: prolactin gene (old: mating preference (female trait))
  u64b     *f1;                           // d; gene-drive: prolactin gene (old: mating preference (female trait))
  double    f;                            // Scaled prolactin (mating preference (female trait)
  u64b     *z0;                           // d; sex chromosomes (old: neutral characters)
  u64b     *z1;                           // d; sex chromosomes (old: neutral characters)
  char       s;                            // Sex of the individual (gene-drive: based on neutral gene)
  char      age;                          // age of individual
  u64b      id;                           // Individuals id (index into ->Individuals)
  //double    TRoam[_RRANGE * _RRANGE];   // individuals time spent profile in the roaming range (NOT USED)
} Individual, *pIndividual;

typedef struct {
  pIndividual       ind;                    // the other individual
} OthersTime,      *pOthersTime;

// Structure to represent a patch
typedef struct patch {
  pIndividual       i[_MAX_INDIVIDUALS];
  char              lat;
  char              lon;
  u64b              ni;                               // The number of individuals currently in this patch
  u64b              n;                                // This patch's "niche value"
  byte              l[_k];                            // The patch traits
  //u64b            nne;                              // Number of neighboring patches (useful, used for ClusterIt())
  //struct patch   *ne[9];                            // Array of neighboring patches (+1 for loop structure)
  u64b              ddne;                             // no of dispersal patches
  struct patch     *dne[(_DRANGE * _DRANGE) + 1];     // array of dispersal patches (plus 1 as above)
#if _DISTANCEDENSITYDISPERSAL
  u64b              disdist[(_DRANGE * _DRANGE) + 1]; // distance data
  double            density;                          // density of this patch
#endif
  OthersTime        Other[_DRANGE * _DRANGE * _MAX_INDIVIDUALS];
  u32b              nothers;                          // no of others spending part of their time in the focal patch
  //char            included;                         // whether it is included in the current patch or not
  //int             rangeid;
  u64b              SpeciesID;                        // Patch's "species value" like "niche value"
  char              AllSpecies[8];                    // the species ids present in the patch
  int               speciespops[8];                   // pop size of each species in the patch
  //int             SRangeIncluded[16];               // for range calculation
  //int             edge;                             // (boolean) whether the patch is at the edge or a center
  } Patch, *pPatch;

// Define a type for our "space"
typedef Patch **Space, ***pSpace;

// This bundles all the data a thread needs into a nice little structure
typedef struct {
  Space           Current,Next;       // Current and next generation spaces
  Individual     *ICurrentData;       // Enough storeage for MAXPOP individuals
  Individual     *INextData;          // Enough storeage for MAXPOP individuals
  pIndividual    *ICurrent;           // Array of pointers to individual storeage
  pIndividual    *INext;              // Array of pointers to individual storeage
  pIndividual    *ICurrentCleared;    // "Cleared" or initial allocation
  pIndividual    *INextCleared;       // "Cleared" or initial allocation
  u64b            nICurrent;          // Number of individual structures in use
  u64b            nINext;             // Number of individual structures in use
  int             id;                 // Thread id
  distribution   *CrossDist[6];       // Byte-wise crossover distribution
  //distribution *MutateDist[6][256]; // Mutation distribution for each character
  //distribution *MutateSkip[6];      // Skip distribution for each character
} ThreadData, *pThreadData;

#endif
