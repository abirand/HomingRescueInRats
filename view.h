#include "types.h"
#ifndef _VIEW_H
#define _VIEW_H


#define GENEDRIVEMODEL        1


// Strategies
#define HOMINGDRIVE           0
#define HOMER                 0
#define HOMER_HS              0
#define HOMER_KO              1

// ********
// important! change these as well
#define INOC_INDIVIDUALS      8

// number of generations the gene-drive individuals get inoculated;
// adjust according to the entry _InocGen below
#define NoInoc                1

// gene-drive individuals inoculation generations
#define FirstInocGen          14

// number of patches the gene-drive individuals get inoculated; adjust according to
#define InocPatchNo           16

#define MATINGCYCLES          6

// ********

// Scale deme graph relative to the maximum deme size from beginning
#define SCALING               0

// Set this to 1 to draw only the basic graphs
#define SIMPLE                0

// Set this to 1 to write pretty post-script files instead of XWD dumps
#define PRETTY_PS             1

// Set this to 1 for gray scale resources
#define BLACK_WHITE           1

// Minimum density for a marker bucket to be considered further for counting as a species
#define SPECIES_K_MINIMUM     .01

// Minimum percent of maximum population needed in a niche for it to be considered "filled"
#define NICHE_FILLED_MINIMUM  .05

// Set this to 1 to use simple colors (colors 1-4)
#define PRESENTATION_MODE     0

// to include distance data
#define READ_DISTANCE         1

// Color modes
#define COLOR_PREF            0
#define COLOR_MARKER          1

// (boolean) this is to check for new cramer's v and species vs preference matrix
#define _NEWRM                1

#ifndef __64b
#define __64b

// These should be reasonably portable 64-bit integer types
//typedef unsigned long long int u64b;
//typedef signed   long long int s64b;
#endif

typedef unsigned long Ulong;
typedef unsigned char Byte;

//#define VERBOSE  // make window creation and file loading verbose
#define PI ((double)3.1415926535897932384)

// Number of possible niches
#define NICHES nTraits

typedef struct {
  Display             *dpy;
  int                  screen;
  Window               win;
  GLXContext           ctx;
  XSetWindowAttributes attr;
  int                  x, y;
  unsigned int         width, height;
  unsigned int         depth;
  GLuint               font;
  int                  id;            // To help the user track the window
} GLWindow;

// Attribute lists for single and double buffered windows
static int attrListSgl[] = {
  GLX_RGBA, GLX_RED_SIZE, 4,
  GLX_GREEN_SIZE, 4,
  GLX_BLUE_SIZE, 4,
  GLX_DEPTH_SIZE, 16,
  None
};
static int attrListDbl[] = {
  GLX_RGBA, GLX_DOUBLEBUFFER,
  GLX_RED_SIZE, 4,
  GLX_GREEN_SIZE, 4,
  GLX_BLUE_SIZE, 4,
  GLX_DEPTH_SIZE, 16,
  None
};

// Simulation output
typedef struct {
  u64b        ni;        // The number of individuals in this patch
  u64b        pnf;       // The niche prefered by most individuals of this patch (fitness)
  double      acp;       // Average carying capactity
  double      apr;       // Average preference
  double      amf;       // Average degree of matching with current patch (fitness)
  double      amm;       // Average degree of matching with current patch (migration)
  double     *af;        // Average traits (fitness)
  double     *am;        // Average traits (migration)
  double     *vf;        // Variance of traits (fitness)
  double     *vm;        // Variance of traits (migration)
  double      m;         // Average mating (choosiness in females)
  double      k;         // Average marker (in males)
  double      f;         // Average female preference
  u64b        n;         // This patch's "niche value"
  u64b        sid;       // This patch's "species value"
  Byte       *l;         // Binary representation of this patch's "niche value" ->n
} patch;

typedef struct {
  int        *idspecies; // species id
  int        *noisolate; // number of ranges (isolated ranges)
  int        *dranges;   // diversity of range sizes
  int        *mrange;    // max range size
  //int      *mdura;     // max duration of species
  double     *meanr;     // mean range size
  double     *sdevr;     // standard deviation
  double     *sir;       // simpson's evenness index
  double     *shr;       // shannon's evenness index
  int        *sppop;     // pop size of species
} rangest;

typedef struct {
  int        *sphno;     // species id
  int        *sphid;     // number of ranges (isolated ranges)
  int        *sphsize;   // diversity of range sizes
  int         msphsize;
} sphistst;

typedef struct {
  int        *idpatch;      // niche type
  int        *noisolpatch;  // number of ranges (isolated ranges)
  int        *dpatchsize;   // diversity of range sizes
  int        *mpatchs;      // max range size
  double     *meanpatchs;   // mean range size
  double     *sdevp;        // standard deviation
  double     *sip;          // simpson's evenness index
  double     *shp;          // shannon's evenness index
  double      am;           // average of max'es
  // the following are for the entire landscape rather than per resource
  int         check;        // niche type
  int         tcount;       // number of ranges (isolated ranges)
  int         maxp;         // diversity of range sizes
  double      meanp;        // max range size
  double      tsd;          // mean range size
  double      tsip;         // simpson's evenness index
  double      tshp;         // shannon's evenness index
  double      frag;         // % patches that are edges
} habitatst;

typedef struct {
  // Properties
  u64b        g;         // Actual generation number (b.c. SKIP)
  u64b        Popsize;   // Total population size as reported by nIndividuals()
  // Per patch stats
  patch     **p;         // Patch space
  // Array / Matrix based stats
  u64b      **Prefs;     // Preference count graph data
  u64b        mPref;     // Value of largest enrty in Prefs
  u64b       *Nichef;    // Niche frequency array
  u64b        mNichef;   // Value of largest entry in Nichef
#if _SPVY
  u64b      **spvy;      // New preference count graph data
  //u64b      mspvy;     // Value of largest enrty in Prefs
#endif
#if _NEWRM
  u64b      **newrm;     // New preference count graph data
  u64b        mnewrm;    // Value of largest enrty in Prefs
#endif
  u64b      **mvk;       // Mating vs. Marker matrix
  u64b        mmvk;      // Max value in mvk matrix
  u64b      **kpn;       // Marker per Niche matrix
  u64b        mkpn;      // Max value in kpn matrix
  u64b      **fpn;       // Female-marker per Niche matrix
  u64b        mfpn;      // Max value in fpn matrix
  u64b      **spn;       // Marker per species matrix
  u64b        mspn;      // Max value in spn matrix
  u64b      **sfpn;      // female-marker per species matrix
  u64b        msfpn;     // Max value in sfpn matrix

  // Species stats per niche
  int        *ns;        // Number of species for each niche
  int         tns;       // Total number of species
  int         fn;        // Number of filled niches
  char       *filled;    // Filled flag for each niche
  double      spfn;      // Number of species per filled niche
  double      vark;      // Average variance of the even entries in the KPN colums

  // Species stats per ecological species
  int        *nss;       // Number of sexual species per ecological species
  int         tnss;      // Total number of sexual species species
  double      sexsps;    // Number of sexual species

  // Stats
  u64b        rMax;      // Largest per-patch population size
  u64b        gMax;      // Largest genotype (population) size
  

  // modified for TALLELE option
  u64b        geno_ze;   // pop sizes for all genotypes
  u64b        geno_on;   // pop sizes for all genotypes
  u64b        geno_tw;   // pop sizes for all genotypes
  u64b        geno_th;   // pop sizes for all genotypes
  u64b        geno_fo;   // pop sizes for all genotypes
  u64b        geno_fi;   // pop sizes for all genotypes
  u64b        geno_si;   // pop sizes for all genotypes
  u64b        geno_se;   // pop sizes for all genotypes
  // original version
  //double    cV;        // Cramer's V
  //u64b      rs;        // Number of "real" species
  //double    pf;        // Average matching (prefered fitness)
  //double    pm;        // Average matching (prefered migration)
  //double    cf;        // Average matching (current fitness)
  //double    cm;        // Average matching (current migration)

  //double    apr;       //
  double      acp;       // Average prolactin 'f'
  double      amc;       // Average gene-drive 'm'
  double      amk;       // Average cat 'k'
  int         off;       // Offspring numbers
  int         mor;       // mortality
  int         spat;      // single paternity
  int         mpat;      // multiple paternity
  
  // Matching / preferences
  u64b        Mpnf;      // The most prefered niche (fitness)
  u64b        Mpnm;      // The most prefered niche (migration)
  
  // range info
  rangest     range;
  habitatst   hab;
  sphistst    sphist;

  // distance info
#if READ_DISTANCE
  u64b      **dvd;       // distance vs. density matrix
  u64b        mdvd;      // Max value in dvd matrix
#endif

  // species coexistence info
  int         spcount;   // no of specialists
  int         gencount;  // no of generalists
  int         tcount;    // total no of species
  int         coexist;   // whether the specialists and generalists coexist (boolean)
  int         multisp;   // whether there are more than one species at one time (boolean)

} Generation;

#endif
