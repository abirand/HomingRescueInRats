#ifndef __PARAMETERS_H
#define __PARAMETERS_H

// this is when running in arbiter cluster (Aaron's)
#define _ARBITER                  0

// this is when running in arbiter cluster (Aaron's)
#define _GENERATIONS              500

// Number of available processors (threads to create); _WIDTH _must_ be divisible by _CPUS
// Note that this version runs only with 1 CPU
#define _CPUS                     1

// Seed for the random number generator, 0 -> seed based on time
#define _RAND_SEED                0

//****************************************
// maximum age of an individual (think in mating cycles!) if generations are overlapping
// (mouse: 12 (2); cats: 10 (5); foxes: 6 (3); rabbits: 12 (3); black rat: 12(2))
#define _MAXAGE                   12

// number of mating cycles/year if generations are overlapping
// has no implications whatsoever, to convert old gens to years
// ***change it in view.h as well!
// (mouse: 6; cats: 2; foxes: 2; rabbits: 4; black rat: 6)
#define _MATINGCYCLES             6

// Birth rate
// (mouse: 6; cats: 5; foxes: 3; rabbits: 4; black rat: 4)
#define _b                        4

//****************************************
// Log() only logs to graph file when (Gen%_SYS) == 0. _SYS == 0 -> disable
#define _SYS                      1

// Log() only logs to graph file when (Gen%_DEME) == 0. _DEME == 0 -> disable
#define _DEME                     1

// Log() only logs to graph file when (Gen%_DISTANCE) == 0. _DISTANCE == 0 -> disable
#define _DISTANCE                 1

// Log() only logs to state file when (Gen%_IND) == 0.  _IND == 0  -> disable
#define _IND                      0
#define _IND_SAMPLE               10

// Width and height of environment
// DO NOT FORGET TO ADJUST _MAXLEVEL ABOVE (32->5; 64->6)
#define _WIDTH                    64
#define _HEIGHT                   64

// Number of different characteristics (for _*e and _*p only)
#define _k                        4

// Maximum number of individuals per patch: this has to do with memory allocation.
// 128 hits the upper max regardless of _KO (even when it is 40, so make it 256!!!)
#define _MAX_INDIVIDUALS          256

// Initial number of individuals to initialize per patch
// watch out for _MAX_INDIVIDUALS for memory and _K0 for dynamics (10 for cats on KI)
#define _INIT_INDIVIDUALS         40

// this is the 'roaming' range for resource search, should be an ODD number
// for larger _DRANGE's RRANGE should be large as well to prevent initial extinctions (5 is OK)
// (13 for cats on KI)
#define _RRANGE                   17

// this is the 'dispersal range' of the offspring.
// should be an ODD number
// previous model used 3, 9, and 15
#define _DRANGE                   17

#define _PANMICTIC                0

//**************************************************************************
#define _FEMALEINOCULATION        1

// default: pC(f) = 0.83
#define _FEMALECUTTINGPROB        0.90

// default: pN(f) = 0.2 (pHOMING = 1-pN=0.8)
#define _FEMALENHEJPROB           0.20

// default: pL(f) = 1.0
#define _FEMALELOFPROB            0.999

// knockout probability of a distal gene, default: pK(f) = 0.83
#define _FEMALEKNOCKOUT           0.77

//------------------------------------------------------------------
#define _MALEINOCULATION          0

// default: pC(m) = 0.29
#define _MALECUTTINGPROB          0.29

// default: pN(m) = 0.97 (pHOMING = 1-pN=0)
#define _MALENHEJPROB             0.97

// default: pL(m) = 1.0
#define _MALELOFPROB              0.999

// knockout probability of a distal gene, default: pK(m) = 1.0
#define _MALEKNOCKOUT             0.29

//**************************************************************************
// Turn approach ONE or TWO 'ON' (Boolean)
//--------------------------------------------------------------------------
// The drive is inserted in haplosufficient essential gene (cis-target),
// generating a loss-of-function mutation with functions:
// if A:    homing drive within haplosufficient female fertility gene
// if B:    homing drive within haplosufficient   male fertility gene
// if C:    A & B together,  then it is a viability locus (no compensation though)
#define _HOMING_A                 0
#define _HOMING_B                 0
#define _HOMING_C                 0
//--------------------------------------------------------------------------
// The drive is inserted in haploinsufficient viability gene, generating a loss-of-function
// mutation but has rescue also targets a distal gene (trans-target) which are:
// if A   (female (AND male) knockout): target is (haplosufficient)     FEMALE fertility
// if B   (female (AND male) knockout): target is (haplosufficient)       MALE fertility
// if A&B (female (AND male) knockout): target is (haplosufficient) 'viability'
// knockout could be in both sexes (adjust _MALEKNOCKOUT above)
#define _HOMERKO_A                1
#define _HOMERKO_B                0
#define _HOMERKO_C                0
// _HOMERKO_x_EXP has two linked copies of the drive which makes MALES with 2 copies sterile
// in this model framework, basically and male with allele '2' are sterile
// because they are tightly linked
//**************************************************************************
// male sterility observed in experimental results when they carry (two) copies of the drive (on the same chromosome)
#define _EXP_MALE_STERILITY       0
//**************************************************************************

// Boolean() only inoculates drive carrying individuals in THREAD() in thread.c if ON
#define _INOC                     1

// the following five should be re-defined in view.h for them to be plotted in view.c
// Initial number of individuals to inoculate per patch
// watch out for _MAX_INDIVIDUALS + _K0 for and dynamics
// 4 for cats on KI
#define _INOC_INDIVIDUALS         1

#define _InocPatchNo              256
#define _InocCoor {2,2,6,2,10,2,14,2,18,2,22,2,26,2,30,2,34,2,38,2,42,2,46,2,50,2,54,2,58,2,62,2,2,6,6,6,10,6,14,6,18,6,22,6,26,6,30,6,34,6,38,6,42,6,46,6,50,6,54,6,58,6,62,6,2,10,6,10,10,10,14,10,18,10,22,10,26,10,30,10,34,10,38,10,42,10,46,10,50,10,54,10,58,10,62,10,2,14,6,14,10,14,14,14,18,14,22,14,26,14,30,14,34,14,38,14,42,14,46,14,50,14,54,14,58,14,62,14,2,18,6,18,10,18,14,18,18,18,22,18,26,18,30,18,34,18,38,18,42,18,46,18,50,18,54,18,58,18,62,18,2,22,6,22,10,22,14,22,18,22,22,22,26,22,30,22,34,22,38,22,42,22,46,22,50,22,54,22,58,22,62,22,2,26,6,26,10,26,14,26,18,26,22,26,26,26,30,26,34,26,38,26,42,26,46,26,50,26,54,26,58,26,62,26,2,30,6,30,10,30,14,30,18,30,22,30,26,30,30,30,34,30,38,30,42,30,46,30,50,30,54,30,58,30,62,30,2,34,6,34,10,34,14,34,18,34,22,34,26,34,30,34,34,34,38,34,42,34,46,34,50,34,54,34,58,34,62,34,2,38,6,38,10,38,14,38,18,38,22,38,26,38,30,38,34,38,38,38,42,38,46,38,50,38,54,38,58,38,62,38,2,42,6,42,10,42,14,42,18,42,22,42,26,42,30,42,34,42,38,42,42,42,46,42,50,42,54,42,58,42,62,42,2,46,6,46,10,46,14,46,18,46,22,46,26,46,30,46,34,46,38,46,42,46,46,46,50,46,54,46,58,46,62,46,2,50,6,50,10,50,14,50,18,50,22,50,26,50,30,50,34,50,38,50,42,50,46,50,50,50,54,50,58,50,62,50,2,54,6,54,10,54,14,54,18,54,22,54,26,54,30,54,34,54,38,54,42,54,46,54,50,54,54,54,58,54,62,54,2,58,6,58,10,58,14,58,18,58,22,58,26,58,30,58,34,58,38,58,42,58,46,58,50,58,54,58,58,58,62,58,2,62,6,62,10,62,14,62,18,62,22,62,26,62,30,62,34,62,38,62,42,62,46,62,50,62,54,62,58,62,62,62}
 

// number of generations the gene-drive individuals get inoculated;
// adjust according to the entry _InocGen below
#define _NoInoc                   1

// gene-drive individuals inoculation generations
#define _InocGen                  {12}

// **********
// (Boolean) NEW flag for converting _Ln to sex chromosomes without deleting the old code
// should be ON
#define _SEXCHROMOSOME            1

//**************************************************************************
// (cats: 0.65; rabbits: 0.7; )
// (KI survival cats: 0.65, then K=3)
#define _SURVIVALPROBABILITY      0.62

// Maximum carrying capacity
// (mice: #K=round(-50*x + 45.75))
// (rats: #K=round(-35*x + 33.78))
#define _K0                       12

// (Boolean) set one of them to 1, rest to zero
// _RAND_DISPERSAL: random dispersal
// _PREF_DISPERSAL: preference based dispersal (based on trait 'y') (not used in gene-drive model)
// _DISTANCEDENSITYDISPERSAL: density dependent dispersal based on distance, define density threshold as well
#define _RAND_DISPERSAL           0
#define _PREF_DISPERSAL           0
#define _DISTANCEDENSITYDISPERSAL 1
#define _DISPERSALCOEF            1
#define _DENSITYCOEF              1

//-------------------------------------------------------------------------------------
// prob. of BEHAVIORAL POLYANDRY per female. This is also the actual GENETIC POLYANDRY, since for simplicity
// BEHAVIORAL POLYANDRY = GENETIC POLYANDRY, ALL multiple matings results in sires "multiple paternity opportunity"
// multiple mating = two males only
#define _POLYANDRYPROB            0.68

// in multiple mating, the advantage of the first male's sperm (early sperm advantage, first come first serve, copulatory plug)
// if both the males have the same genotype (w,w or dr,dr), then m1 = 0.7, m2 = 1-m1 = 0.3
#define _FIRSTSPERMADVANTAGE      0.5

// this is the 'offset' if the males are (w,dr), if (w,dr)=(0.7+0.2, (1-0.7)-0.2); if if (dr,w)=(0.7-0.2, (1-0.7)+0.2);
// THE TWO ABOVE CANNOT BE GREATER THAN 1
#define _SPERMCOMPETITIONCOEF     0.0

//-------------------------------------------------------------------------------------

// Number of loci per characteristic; _must_ be even if _RANDOM_GENOTYPE == 0
// the warning above about being even is for the crossover to work,
// _Le = 1 works (in terms of the code, no problem, ok if you don't want crossover)
// _Ln = 1 in order to be used as _SEXCHROMOSOME
// _Lm = 8 used as a drive (8 is necessary in GeneDrive() in thread.c)
// _Lf = 8 used as a prolactin (8 may not be necessary in GeneDrive() in thread.c)
#define _Le                       1
#define _Lp                       1
#define _Lm                       1
#define _Lk                       1
#define _Lf                       1
#define _Ln                       1

// Number of bits per allele; _must_ be a member of {1,2,4,8}
// usually all are set to 1, and _An is set to 8
#define _Ae                       1
#define _Ap                       1
#define _Am                       1
#define _Ak                       1
#define _Af                       4
#define _An                       1

// Mutation rate (usually all are set to 0.00001, _mun set to 0.001)
// cannot be zero
#define _mue                      0.00001
#define _mup                      0.00001
#define _mum                      0.00001
#define _muk                      0.00001
#define _muf                      0.00001
#define _mun                      0.00001

// You _must_ set these to 0 turn off mutation
#define _use_mue                  0
#define _use_mup                  0
#define _use_mum                  0
#define _use_muk                  0
#define _use_muf                  0
#define _use_mun                  0

// When enabled, this will cause niche to use the faster "skip" method to do mutation
#define  SKIP_MUTATION            0

// Crossover rate
// this is corrected on 02.07.20 to include 6 entries, it was erroronously 3 in the original
// model (6 for 6 traits: x,y,m,f,k,n)
// these are probabilities, not r=0.5. could be set to different values
// original code treats loci independently, there is free recombination between the loci/trait
// in gene-drive model, I turned the cross over off, traits are independent though.
// could be turned on, but NOT in drive related loci
//#define _Xi                     {.5,.5,.5,.5,.5,.5}
#define _Xi                       {0.,0.,0.,0.,0.,0.}
//-------------------------------------------------------------------------------------

// When enabled, this will cause niche to use Micheal's faster version of exp()
#define  _FAST_EXP                1

// Will cause some extra printing to be done if == 1 (2 for more)
#define _VERBOSE                  1

// Causes some timing information to be computed and printed
#define _BENCHMARK                1

// Prefix to use for output file names (before .ind, etc).  See FormatName() for more details
#define _FILE_NAME                "%g"

#endif
