#include <stdio.h>
#include <gsl/gsl_fit.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/select.h>
#include <pthread.h>
#include <unistd.h>
#include <GL/glx.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "view.h"
#include "types.h"
#if PRETTY_PS
#include "ps.h"
#endif

/***************************************************************************************
 * Globals
 ***************************************************************************************/

Generation      *Gen;                // pointers to generations
u64b             nGen;               // number of generations
u64b             Width,Height;       // Width and Height of the patch space
int              nTraits;            // Number of traits (_k)
double           Max;                // _K0
double           _b;
int              GlobalMax;          // global max population size to scale and record
int              GlobalMin;          // global min population size to record
int              GenMax;             // overall genotype pop max to scale
int              PPopMax;            // overall patch pop max to scale
int              Clean;              // True if the data file contains a trailer
int              Amcg=-1;            // The generation which fisrt saw amc >= .5
int              Amcf=-1;            // The frame which fisrt saw amc >= .5

char             *FNPrefix;

int              Sexual;             // Flags for the sexual options
int              Selection;          // Flags for the sexual options
int              Speciesid;          // Flags for the sexual options
#if READ_DISTANCE
int              DVDy;               // Max distance, y dimension of dvd matrix
#endif
float            Trans[2],Scale=22;  // graph translation and scale
float            HackScale;          // graph translation and scale
int              Mx=10000,My=10000;  // Mouse coordinates
int              Speed = 10, Frame;  // FPSm and currently drawn generation
int              X=1000,Y=900;       // Width and height of window
int              dump;               // dump mode flag
int              ColorMode;          // Draw colors for pref or colors for marker
int              Filter[1024];       // Filters which niches are drawn

int              Timetocolonize = 0;
int              TPopSizeReduction = 0;  // time when popsize falls below 1% of maximum popsize
int              TMaxPop = 0;            // time when max popsize is reached
int              TMinPop = 0;            // time when max popsize is reached

// Simple vector type
typedef struct {
  float          x;
  float          y;
  float          z;
} vec3f_t;

typedef struct {
  s32b           ox;        // Mouse down position
  s32b           oy;
  s32b           nx;        // Current mouse position
  s32b           ny;
  u32b           d;         // Mouse down flag
  vec3f_t        rot;       // Rotation vector
  vec3f_t        trn;       // Translation vector
  s32b           x;
  s32b           y;
  s32b           w;
  s32b           h;
} plot_3d_t;

plot_3d_t        Plot;


#if PRETTY_PS
FILE            *PSFile0;
FILE            *PSFile1;
#endif

/* Color blend array */
Byte RB[78][3] = {{0xff, 0x00, 0x00},{0xff, 0x10, 0x00},{0xff, 0x1f, 0x00},{0xff, 0x2f, 0x00},
		  {0xff, 0x3e, 0x00},{0xff, 0x4e, 0x00},{0xff, 0x5d, 0x00},{0xff, 0x6d, 0x00},
		  {0xff, 0x7c, 0x00},{0xff, 0x8c, 0x00},{0xff, 0x9c, 0x00},{0xff, 0xab, 0x00},
		  {0xff, 0xbb, 0x00},{0xff, 0xca, 0x00},{0xff, 0xda, 0x00},{0xff, 0xe9, 0x00},
		  {0xff, 0xf9, 0x00},{0xf6, 0xff, 0x00},{0xe6, 0xff, 0x00},{0xd6, 0xff, 0x00},
		  {0xc7, 0xff, 0x00},{0xb7, 0xff, 0x00},{0xa8, 0xff, 0x00},{0x98, 0xff, 0x00},
		  {0x89, 0xff, 0x00},{0x79, 0xff, 0x00},{0x6a, 0xff, 0x00},{0x5a, 0xff, 0x00},
		  {0x4a, 0xff, 0x00},{0x3b, 0xff, 0x00},{0x2b, 0xff, 0x00},{0x1c, 0xff, 0x00},
		  {0x0c, 0xff, 0x00},{0x00, 0xff, 0x03},{0x00, 0xff, 0x13},{0x00, 0xff, 0x22},
		  {0x00, 0xff, 0x32},{0x00, 0xff, 0x42},{0x00, 0xff, 0x51},{0x00, 0xff, 0x61},
		  {0x00, 0xff, 0x70},{0x00, 0xff, 0x80},{0x00, 0xff, 0x8f},{0x00, 0xff, 0x9f},
		  {0x00, 0xff, 0xae},{0x00, 0xff, 0xbe},{0x00, 0xff, 0xce},{0x00, 0xff, 0xdd},
		  {0x00, 0xff, 0xed},{0x00, 0xff, 0xfc},{0x00, 0xf2, 0xff},{0x00, 0xe3, 0xff},
		  {0x00, 0xd3, 0xff},{0x00, 0xc4, 0xff},{0x00, 0xb4, 0xff},{0x00, 0xa4, 0xff},
		  {0x00, 0x95, 0xff},{0x00, 0x85, 0xff},{0x00, 0x76, 0xff},{0x00, 0x66, 0xff},
		  {0x00, 0x57, 0xff},{0x00, 0x47, 0xff},{0x00, 0x38, 0xff},{0x00, 0x28, 0xff},
		  {0x00, 0x18, 0xff},{0x00, 0x09, 0xff},{0x07, 0x00, 0xff},{0x16, 0x00, 0xff},
		  {0x26, 0x00, 0xff},{0x35, 0x00, 0xff},{0x45, 0x00, 0xff},{0x54, 0x00, 0xff},
		  {0x64, 0x00, 0xff},{0x74, 0x00, 0xff},{0x83, 0x00, 0xff},{0x93, 0x00, 0xff},
		  {0xa2, 0x00, 0xff},{0xb2, 0x00, 0xff}};

#if PRESENTATION_MODE
/* simple color array */
Byte simple[8][3] = {{0,255,0},{255,0,0},{0,0,255},{255,255,255},{255,0,255},{255,255,0},
		     {0,255,255},{255,110,0}};
#endif

void Black()  { glColor3ub(  0,   0,   0); }
void Cyan()   { glColor3ub(  0, 255, 255); }
void Yellow() { glColor3ub(255, 255,   0); }
void Red()    { glColor3ub(255,   0,   0); }
void Purple() { glColor3ub(255,   0, 255); }
void Blue()   { glColor3ub(  0,   0, 255); }
void White()  { glColor3ub(255, 255, 255); }
void Green()  { glColor3ub(  0, 255,   0); }


/***************************************************************************************
 * Font / GL / X / etc
 ***************************************************************************************/

void ViewPort2D(GLWindow *glw)
{
  glViewport(0, 0, glw->width, glw->height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, glw->width, glw->height, 0.0, 1, -100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void ViewPort3D(int x, int y, int w, int h)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(x, y, w, h); 
  gluPerspective(55.0f,(float)w/(float)h,0.8f,10.0f);
  glTranslatef(0.0f, 0.0f, -2.0f);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

static GLuint BuildFont(GLWindow *glw, char *fname)
{
  XFontStruct *font;  GLuint fbase;   

  /* Storage for 96 characters */
  fbase = glGenLists(96);

  /* Load a font with a specific name in "Host Portable Character Encoding" */
  if ( !(font = XLoadQueryFont(glw->dpy, fname)) ) {
    fprintf(stderr,"Could not load font \"%s\"\n",fname);
    exit(-1);
  }

  /* Build 96 display lists out of our font starting at char 32 */
  glXUseXFont(font->fid, 32, 96, fbase);
  XFreeFont(glw->dpy, font);

  return fbase;
}

/*
static void DeleteFont(GLuint font)
{
  glDeleteLists(font, 96);
}
*/

static void printGLf(GLuint font, const char *fmt, ...)
{
  va_list ap;  char text[256];

  if (fmt == NULL) return;

  va_start(ap, fmt);  
  vsprintf(text, fmt, ap);
  va_end(ap);
  glPushAttrib(GL_LIST_BIT);
  glListBase(font - 32);
  glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
  glPopAttrib();
}

/* called when window is resized */
static void ResizeGLScene(unsigned int width, unsigned int height)
{
  glViewport(0, 0, width, height); 
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0., width, height, 0., 1, -100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

/* destroy window */
static void KillGLWindow(GLWindow *glw)
{
  if (glw->ctx) {
    glDeleteLists(0, 65535);
    if (!glXMakeCurrent(glw->dpy, None, NULL))
      printf("Could not release drawing context.\n");
    glXDestroyContext(glw->dpy, glw->ctx);
    glw->ctx = NULL;
  }
  XCloseDisplay(glw->dpy);
}

/* create window */
static void CreateGLWindow(GLWindow *glw, char* title, int width, int height)
{
  XVisualInfo *vi;  Colormap cmap;  Window twin;  Atom wmDelete;  int glxM, glxm, t;

  /* Connect to X */
  if( !(glw->dpy    = XOpenDisplay(0)) ) {
    fprintf(stderr,"Cannot connect to X server!\n");
    exit(1);
  }
  glw->screen = DefaultScreen(glw->dpy);

  /* Get the appropriate visual */
  glXQueryVersion(glw->dpy, &glxM, &glxm);
#ifdef VERBOSE
  printf("%d: glX-Version %d.%d\n", glw->id, glxM, glxm);
#endif
  if ( (vi = glXChooseVisual(glw->dpy, glw->screen, attrListDbl)) )
#ifdef VERBOSE
    printf("%d: Selected doublebuffered mode.\n", glw->id);
#else
  ;
#endif
  else {
    vi = glXChooseVisual(glw->dpy, glw->screen, attrListSgl);
#ifdef VERBOSE
    printf("%d: Selected singlebuffered mode.\n", glw->id);
#endif
  }

  /* Create a GLX context and color map */
  glw->ctx               = glXCreateContext(glw->dpy, vi, 0, GL_TRUE);
  cmap                   = XCreateColormap(glw->dpy, RootWindow(glw->dpy, vi->screen), vi->visual, AllocNone);
  glw->attr.colormap     = cmap;
  glw->attr.border_pixel = 0;

  /* Create the window */
  glw->attr.event_mask = ExposureMask | KeyPressMask | ButtonPressMask| ButtonReleaseMask | StructureNotifyMask | PointerMotionMask;
  glw->win = XCreateWindow(glw->dpy, RootWindow(glw->dpy, vi->screen),
                           0, 0, width, height, 0, vi->depth, InputOutput, vi->visual,
                           CWBorderPixel | CWColormap | CWEventMask, &glw->attr);
  wmDelete = XInternAtom(glw->dpy, "WM_DELETE_WINDOW", True);
  XSetWMProtocols(glw->dpy, glw->win, &wmDelete, 1);
  XSetStandardProperties(glw->dpy, glw->win, title, title, None, NULL, 0, NULL);
  XMapRaised(glw->dpy, glw->win);

  /* Connect the glx-context to the window */
  glXMakeCurrent(glw->dpy, glw->win, glw->ctx);
  XGetGeometry(glw->dpy, glw->win, &twin, &glw->x, &glw->y,
              &glw->width, &glw->height, (unsigned int*)&t, &glw->depth);
#ifdef VERBOSE
  printf("%d: Depth %d\n", glw->id, glw->depth);
  if (glXIsDirect(glw->dpy, glw->ctx)) printf("%d: Direct Rendering enabled.\n", glw->id);
  else                                 printf("%d: Direct Rendering disabled.\n", glw->id);
#endif
}

void yeild()
{
  struct timeval tv;

  // Sleep for a ms
  tv.tv_sec  = 0;
  tv.tv_usec = 30000;
  select(0, NULL,  NULL, NULL, &tv);
}

static int AnyEvent(Display *d, XEvent *e, XPointer p)
{
  return True;
}

/***************************************************************************************
 * Data file loading code
 ***************************************************************************************/

#define ALLOC(x,y)  if (!(x = malloc(y))) { perror("niche_view"); exit(2); }
#define ZALLOC(x,y) if (!(x = malloc(y))) { perror("niche_view"); exit(2); } else memset(x,0,y);

int ReadSys(char *fn)
{
  u64b i,j;  FILE *f;  char b[1024];  int lm,g=0;  double min;
  GlobalMax = 0;
  GlobalMin = 200000;
  GenMax = 0;
  PPopMax = 0;
  // Open file and read the header
  Clean = 1;
  if (!(f = fopen(fn,"r"))) {
    sprintf(b,"niche_view: .sys-open: \"%s\"",fn);
    perror(b);
    return 1;
  }
  fseek(f,i=-1,SEEK_END);
  while( (lm=fgetc(f)) != '\n' ) {
    if( lm == EOF ) {
      fprintf(stderr,"niche_view: data file looks VERY wrong..\n");
      return 3;
    }
    fseek(f,--i,SEEK_END);
  }
  if( fscanf(f,"%llu", &nGen) != 1 ) { 
    fprintf(stderr, "niche_view: sys: trailer error\n");
    fprintf(stderr, "niche_view: sys: attempting to guess... (500) \n");
    nGen = 500;
    g = 1;
  }
  rewind(f);
  if (fscanf(f,"%d %lf %lf %llu %llu\n",&nTraits,&Max,&_b,&Width,&Height) != 5) 
    { fprintf(stderr, "niche_view: sys: header error\n"); return 1; }
  if(nTraits > 64) {
    fprintf(stderr, "niche_view: sys: _k too large.  (0 < _k < 64)\n");
    return 1;
  }
  // Auto-Scale
  Scale = ((((X-280)/Width)<((Y-410)/(Height)))?((X-280)/Width):((Y-410)/(Height)));
  HackScale = Scale;
#ifdef VERBOSE
  fprintf(stderr,"Sys Header: %d traits, %u Generations, %llux%llu, %lf Max\n",nTraits,nGen,Width,Height,Max);
#endif
  // Now that the header is read in, compute the minimum number of species
  // needed to call a niche "filled"
  min = NICHE_FILLED_MINIMUM*((Width*Height)/(nTraits))*Max;
  // Allocate generation structures
  ZALLOC(Gen,  nGen*sizeof(Generation));
  for(i=0; i<nGen; i++) {
    // Allocate patch data structures
    ALLOC(Gen[i].p,             Width*sizeof(patch*));
    for(j=0; j<Width; j++)
      ZALLOC(Gen[i].p[j],       Height*sizeof(patch));
    // Allocate preference matrix
    ALLOC(Gen[i].Prefs,         NICHES*sizeof(u64b*));
    for(j=0; j<NICHES; j++)
      ZALLOC(Gen[i].Prefs[j],   NICHES*sizeof(u64b));
    // Allocate niche frequency array
    ZALLOC(Gen[i].Nichef,       NICHES*sizeof(u64b));
    // Allocate species count array
    ZALLOC(Gen[i].ns,           NICHES*sizeof(int));
    // Allocate fille niche flags
    ZALLOC(Gen[i].filled,       NICHES*sizeof(char));
    ZALLOC(Gen[i].nss,          16*sizeof(int));
  }

  // Read system data
  for(i=0; i<nGen; i++) {
    // Read generation header
    if( fscanf(f, "gen: %llu %llu ", &Gen[i].g, &Gen[i].Popsize) != 2 ) {
      fprintf(stderr,"niche_view: sys: generation header error (%llu)\n",i);
      if( g == 1 ) {
        nGen = i;
        Clean = 0;
        fprintf(stderr,"niche_view: sys: trailer was not present; using (%llu) as nGen instead.\n",i); 
        break;
      } else {
        return 1;
      }
    }
    // If empty gen, skip it
    // this skipping option didn't work when extinct?
    //if( !Gen[i].Popsize )
    //  continue;
    
    // Read generation data
    // in gene-drive model acp is changed to i->f, i.e. prolactin

    if( fscanf(f, "%lf %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu %lf %lf %d %d %d %d ",
              &Gen[i].acp, &Gen[i].gMax, &Gen[i].rMax, &Gen[i].geno_ze, &Gen[i].geno_on,
              &Gen[i].geno_tw, &Gen[i].geno_th, &Gen[i].geno_fo, &Gen[i].geno_fi, &Gen[i].geno_si, &Gen[i].geno_se, &Gen[i].amc, &Gen[i].amk, &Gen[i].off, &Gen[i].mor, &Gen[i].mpat, &Gen[i].spat )  != 17 )
    { fprintf(stderr,"niche_view: sys: generation error (%llu)\n",i); return 1; }
    // to scale properly during plotting
    //printf("before if GenMax: %llu\n", GenMax);
    
    if( i <= FirstInocGen ){
      if( Gen[i].Popsize > GlobalMax ){
        GlobalMax = Gen[i].Popsize;
        TMaxPop = i;
        //printf("GlobalMax: %llu, t: %llu\n", GlobalMax, TMaxPop);
      }
    }
    if( i > FirstInocGen ){
      if( Gen[i].Popsize < GlobalMin ){
        GlobalMin = Gen[i].Popsize;
        TMinPop = i;
        //printf("GlobalMax: %llu, t: %llu\n", GlobalMax, TMaxPop);
      }
    }
    if( Gen[i].gMax > GenMax ){
      GenMax = Gen[i].gMax;
      //printf("******if, Gen[i].gMax>GenMax: %llu > %llu\n", Gen[i].gMax, GenMax);
    }
    if( Gen[i].rMax > PPopMax ){
      PPopMax = Gen[i].rMax;
      //printf("******if, Gen[i].gMax>GenMax: %llu > %llu\n", Gen[i].gMax, GenMax);
    }
    //printf("in sys, Gen[i].gMax: %llu, and GenMax: %llu\n\n", Gen[i].gMax, GenMax);
    
    // Read in the terminating newline
    fscanf(f,"\n");
  }
  // Return success
  return 0;
}

int ReadDeme(char *fn)
{
  u64b i,j,k,l,t,g;  FILE *f;  patch *p;  int traits;  u64b width,height;  double max;

  /* Open file and read the header */
  if (!(f = fopen(fn,"r")))
    { perror("niche_view: .deme-open"); return 2; }
  fseek(f,i=-1,SEEK_END);
  while( fgetc(f) != '\n') {
    fseek(f,--i,SEEK_END);
  }
  if( fscanf(f,"%llu", &g) != 1 ) {
    fprintf(stderr, "niche_view: deme: trailer error\n"); 
    fprintf(stderr, "niche_view: deme: using sys value (or guess) instead...\n"); 
    g = 1;
  } else {
    g = 0;
  }
  rewind(f);
  if (fscanf(f,"%d %lf %lf %llu %llu\n",&traits,&max,&_b,&width,&height) != 5) 
    { fprintf(stderr, "niche_view: deme: header error\n"); return 2; }
  if( (traits != nTraits) || (max != Max) || (width != Width) || (height != Height) )
    { fprintf(stderr, "niche_view: deme: deme file does not match sys file (header)\n"); return 2; }
  if(nTraits > 64) {
    fprintf(stderr, "niche_view: deme: _k too large.  (0 < _k < 64)\n");
    return 1;
  }
#ifdef VERBOSE
  fprintf(stderr,"Deme Header: %d traits, %u Generations, %llux%llu, %lf Max\n",nTraits,nGen,Width,Height,Max);
#endif

  /* Read deme data */
  for(i=0; i<nGen; i++) {
    // Read the generation "header" just to make sure no generation 'skew' is introduced.
    if( fscanf(f, "gen: ") == EOF ) {
      fprintf(stderr,"niche_view: deme: generation header error (%llu)\n",i); 
      if( g == 1 ) {
        fprintf(stderr,"niche_view: deme: sys ends at (%llu) but deme ends at (%llu)\n.",nGen,i); 
        fprintf(stderr,"niche_view: deme: setting nGen to (%llu)\n.",i); 
        nGen = i;
      } else {
        return 1;
      }
    }
    for(j=0; j<Width; j++) {
      for(k=0; k<Height; k++) {
        p = &Gen[i].p[j][k];
	/* Allocate needed data structures */
	ALLOC(p->l,      nTraits*sizeof(u64b));
	ALLOC(p->af,     nTraits*sizeof(u64b));
	ALLOC(p->am,     nTraits*sizeof(u64b));
	ALLOC(p->vf,     nTraits*sizeof(u64b));
	ALLOC(p->vm,     nTraits*sizeof(u64b));
	/* Read deme header */
	if( fscanf(f, "%llu %llu ", &l, &p->ni) != 2 )
	  { fprintf(stderr,"niche_view: deme: per-deme generation header error (%llu)\n",i); return 1; }
	/* Read this patch's stats */
	if( fscanf(f, "%llu %lf %lf %lf %lf %llu ", &p->pnf, &p->acp, &p->apr, &p->amf, &p->amm, &p->n) != 6 )
	  { fprintf(stderr,"niche_view: deme: deme data error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	/* Read niche traits */
	for(l=0; l<nTraits; l++) { 
	  if( fscanf(f, "%llu ", &t) != 1)
	    { fprintf(stderr,"niche_view: deme: trait error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	  p->l[l] = ((t)?(1):(0));
	}
	/* Read average traits */
	for(l=0; l<nTraits; l++) 
	  if( fscanf(f, "%lf ", &p->af[l]) != 1 )
	    { fprintf(stderr,"niche_view: deme: average xtrait error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	/* Read variance of traits */
	for(l=0; l<nTraits; l++) 
	  if( fscanf(f, "%lf ", &p->vf[l]) != 1 )
	    { fprintf(stderr,"niche_view: deme: variance xtrait error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	// Read mating and marker if _SEXUAL flag is set
	if( fscanf(f, "%llu ", &t) != 1 )
	  { fprintf(stderr,"niche_view: deme: _SEXUAL flag error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	Sexual = t;
	if( Sexual ) {
	  if( fscanf(f, "%lf ", &p->m) != 1 )
	    { fprintf(stderr,"niche_view: deme: mating error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	  if( fscanf(f, "%lf ", &p->k) != 1 )
	    { fprintf(stderr,"niche_view: deme: mating error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	  // Read female preference if _SEXUAL_SELECTION flag is set
      // while reading deme that is always on for gene-drive this error should never come up
	  if( fscanf(f, "%llu ", &t) != 1 )
	    { fprintf(stderr,"niche_view: deme: _SEXUAL_SELECTION flag error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	  Selection = t;
	  if( Selection ) {
	    if( fscanf(f, "%lf ", &p->f) != 1 )
	      { fprintf(stderr,"niche_view: deme: fpref error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
	  }
        }
        // Read species
        if( fscanf(f, "%llu ", &t) != 1 )
          { fprintf(stderr,"niche_view: deme: species ID flag error (%llu,[%llu,%llu])\n",i,j,k); return 1; }
        if( fscanf(f, "%llu ", &p->sid) != 1 )
        { fprintf(stderr,"niche_view: deme: deme species id and type data error (%llu,[%llu,%llu])\n",i,j,k); return 1; }

        // Read terminating newline
        fscanf(f, "\n");
      }
    }
  }

  // Return success
  return 0;
}

int ReadRange(char *fn)
{
  u64b i,l,g;  FILE *f;  rangest *r;  int traits;  u64b width,height;  double max;

  /* Open file and read the header */
  if (!(f = fopen(fn,"r")))
  { perror("niche_view: .range-open"); return 2; }
  fseek(f,i=-1,SEEK_END);
  while( fgetc(f) != '\n') {
    fseek(f,--i,SEEK_END);
  }
  if( fscanf(f,"%llu", &g) != 1 ) {
    fprintf(stderr, "niche_view: range: trailer error\n"); 
    fprintf(stderr, "niche_view: range: using sys value (or guess) instead...\n"); 
    g = 1;
  } else {
    g = 0;
  }
  rewind(f);
  if (fscanf(f,"%d %lf %lf %llu %llu\n",&traits,&max,&_b,&width,&height) != 5) 
  { fprintf(stderr, "niche_view: range: header error\n"); return 2; }
  if( (traits != nTraits) || (max != Max) || (width != Width) || (height != Height) )
  { fprintf(stderr, "niche_view: range: range file does not match sys file (header)\n"); return 2; }
  if(nTraits > 64) {
    fprintf(stderr, "niche_view: range: _k too large.  (0 < _k < 64)\n");
    return 1;
  }
#ifdef VERBOSE
  fprintf(stderr,"Range Header: %d traits, %u Generations, %llux%llu, %lf Max\n",nTraits,nGen,Width,Height,Max);
#endif

  // Read range data
  for(i=0; i<nGen; i++) {
    // Read the generation "header" just to make sure no generation 'skew' is introduced.
    if( fscanf(f, "gen: ") == EOF ) {
      fprintf(stderr,"niche_view: range: generation header error (%llu)\n",i); 
      if( g == 1 ) {
        fprintf(stderr,"niche_view: range: sys ends at (%llu) but rane ends at (%llu)\n.",nGen,i); 
        fprintf(stderr,"niche_view: range: setting nGen to (%llu)\n.",i); 
        nGen = i;
      } else {
        return 1;
      }
    }

    r = &(Gen[i].range);
    // Allocate needed data structures 
    ALLOC(r->idspecies, ((int)(pow((double)(2),(double)(nTraits))))*sizeof(int));
    ALLOC(r->noisolate, ((int)(pow((double)(2),(double)(nTraits))))*sizeof(int));
    ALLOC(r->dranges,   ((int)(pow((double)(2),(double)(nTraits))))*sizeof(int));
    ALLOC(r->mrange,    ((int)(pow((double)(2),(double)(nTraits))))*sizeof(int));
    ALLOC(r->meanr,     ((int)(pow((double)(2),(double)(nTraits))))*sizeof(u64b));
    ALLOC(r->sdevr,     ((int)(pow((double)(2),(double)(nTraits))))*sizeof(u64b));
    ALLOC(r->sir,       ((int)(pow((double)(2),(double)(nTraits))))*sizeof(u64b));
    ALLOC(r->shr,       ((int)(pow((double)(2),(double)(nTraits))))*sizeof(u64b));
    ALLOC(r->sppop,     ((int)(pow((double)(2),(double)(nTraits))))*sizeof(int));

    u64b g;


    // Read range header
    for(l=0; l<(((int)(pow((double)(2),(double)(nTraits))))+1); l++){
      //printf("\n--------l is %d/%d",l,((int)(pow((double)(2),(double)(nTraits)))));
      fscanf(f, "%llu ", &g);
      //printf("* gen is %lf", g);
      fscanf(f, "%d ",  &r->idspecies[l]);
      //printf("\n id %d ",r->idspecies[l]); 
      fscanf(f, "%d ",  &r->noisolate[l]);
      //printf("noiso %d ", r->noisolate[l]);
      fscanf(f, "%d ",  &r->dranges[l]);
      //printf("dr %d ", r->dranges[l]);
      fscanf(f, "%d ",  &r->mrange[l]);
      //printf("mr %d ", r->mrange[l]);
      fscanf(f, "%lf ", &r->meanr[l]);
      //printf("mean %lf ", r->meanr[l]);
      fscanf(f, "%lf ", &r->sdevr[l]);
      //printf("sd %lf ", r->sdevr[l]);
      fscanf(f, "%lf ", &r->sir[l]);
      //printf("sir %lf ", r->sir[l]);
      fscanf(f, "%lf ", &r->shr[l]);
      //printf("shr %lf ", r->shr[l]);
      fscanf(f, "%d ",  &r->sppop[l]);
      fscanf(f, "\n"); // Read terminating newline
      /*if( fscanf(f, "%d %d %d %d %lf %lf %lf %lf\n",  &r->idspecies[l], &r->noisolate[l],  
          &r->dranges[l], &r->mrange[l], &r->meanr[l], &r->sdevr[l], &r->sir[l], &r->shr[l]) != 8 )
      { fprintf(stderr,"niche_view: range: range error (%llu)\n",i); return 1; }*/
    }
  }

  // Return success
  return 0;
}


#if READ_DISTANCE
int ReadDistance(char *fn)
{
  u64b i,g;  FILE *f;  int j, k, dummygen;

  /* Open file and read the header */
  if (!(f = fopen(fn,"r")))
  { perror("niche_view: .distance-open"); return 2; }
  fseek(f,i=-1,SEEK_END);
  while( fgetc(f) != '\n') {
    fseek(f,--i,SEEK_END);
  }
  // checks for a trailer but doesn't really care what value. hmm...
  if( fscanf(f,"%llu ", &g) != 1 ) {
    fprintf(stderr, "niche_view: distance: trailer error\n");
    fprintf(stderr, "niche_view: distance: using sys value (or guess) instead...\n");
    g = 1;
  } else {
    g = 0;
  }
  rewind(f);
  
  if (fscanf(f,"%d\n",&DVDy) != 1)
  { fprintf(stderr, "niche_view: distance: header error\n"); return 2; }
  
  //printf("nGen: %d\n", nGen);
  /* Read distance data */
  for(i=0; i<nGen; i++) {

    // Read the generation "header" just to make sure no generation 'skew' is introduced.
    if( fscanf(f, "gen: ") == EOF ) {
      fprintf(stderr,"niche_view: distance: generation header error (%llu)\n",i);
      if( g == 1 ) {
        fprintf(stderr,"niche_view: distance: sys ends at (%llu) but distance ends at (%llu)\n.",nGen,i);
        fprintf(stderr,"niche_view: distance: setting nGen to (%llu)\n.",i);
        nGen = i;
      } else {
        return 1;
      }
    }

    Gen[i].mdvd = 0;
    
    if( fscanf(f, "%d ", &dummygen) != 1 )
      { fprintf(stderr,"niche_view: distance: dvd matrix flag error (%llu)\n",i); return 1; }
    
    // Allocate space for this generation's dvd
    ALLOC(Gen[i].dvd, 12*sizeof(u64b));
    for(j=0; j<12; j++)
      ALLOC(Gen[i].dvd[j], DVDy*sizeof(u64b));
    // Read in the mvk matrix
    for(j=0; j<12; j++){
      for(k=0; k<DVDy; k++) {
        // Read value
        fscanf(f, "%llu ", &Gen[i].dvd[j][k]);
        //printf("(%d, %d): %llu ", j, k, Gen[i].dvd[j][k]);
        // Check for max
        if(Gen[i].dvd[j][k] > Gen[i].mdvd)
          Gen[i].mdvd = Gen[i].dvd[j][k];
      }
    // Read terminating newline
    fscanf(f, "\n");
    }
  }

  // Return success
  return 0;
}
#endif


int readDataFile(char *fn)
{
  char sys[1024], deme[1024];  int err;
#if READ_DISTANCE
  char distance[1024];
#endif
  // fn is really just a file prefix, so build 2 real file names
  sprintf(sys,      "%s.sys",     fn);
  sprintf(deme,     "%s.deme",    fn);
  //sprintf(range,    "%s.range",   fn);
#if READ_DISTANCE
  sprintf(distance, "%s.distance",   fn);
#endif
  // Read in system and then deme data
  if( (err=ReadSys(sys))   )          return err;
  //printf("after reading sys, GenMax: %llu\n",GenMax);
  if( (err=ReadDeme(deme)) )          return err;
  //if( (err=ReadRange(range)) )        return err;
#if READ_DISTANCE
  if( (err=ReadDistance(distance)) )  return err;
#endif
  // Return success
  return 0;
}

/***************************************************************************************
 * Util functions used by drawing code
 ***************************************************************************************/

static void InitGLWindow(GLWindow *glw)
{
  glShadeModel(GL_SMOOTH);
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClearDepth(1.0f);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
#if PRESENTATION_MODE
  glw->font = BuildFont(glw, "-adobe-helvetica-bold-r-normal-*-*-180-*-*-p-*-iso10646-1");
#else
  glw->font = BuildFont(glw, "fixed");   
#endif
  ResizeGLScene(glw->width, glw->height);
  glFlush();
}

#if PRESENTATION_MODE
Byte red(double w)
{
  if (w < 0.) w = 0.;
  if (w > 1.) w = 1.;
  return simple[(int)(w*7+.5)][0];
}

Byte green(double w)
{
  if (w < 0.) w = 0.;
  if (w > 1.) w = 1.;
  return simple[(int)(w*7.5)][1];
}

Byte blue(double w)
{
  if (w < 0.) w = 0.;
  if (w > 1.) w = 1.;
  return simple[(int)(w*7.5)][2];
}
#else
Byte red(double w)
{
  if (w < 0.) w = 0.;
  if (w > 1.) w = 1.;
  return RB[(int)(.5+77*(1.-w))][0];
}

Byte green(double w)
{
  if (w < 0.) w = 0.;
  if (w > 1.) w = 1.;
  return RB[(int)(.5+77*(1.-w))][1];
}

Byte blue(double w)
{
  if (w < 0.) w = 0.;
  if (w > 1.) w = 1.;
  return RB[(int)(.5+77*(1.-w))][2];
}
#endif

Byte toByte(double x)
{
  if (x <   0) return 0;
  if (x > 255) return 255;
  return (unsigned char) x;
}

void glOutlineCircle(int x, int y, double r)
{
  int k;

  glBegin(GL_LINE_LOOP);
  for(k=0;k<32;k++)
    glVertex2d(x+r*cos(k*2*PI/32),y+r*sin(k*2*PI/32));
  glEnd();
}

void glCircle(int x, int y, double r)
{
  int k;

  glBegin(GL_POLYGON);
  for(k=0;k<32;k++)
    glVertex2d(x+r*cos(k*2*PI/32),y+r*sin(k*2*PI/32));
  glEnd();
}

/***************************************************************************************
 * Components of the graph
 ***************************************************************************************/

void w()
{
  glColor3ub(255,255,255);
}

void y()
{
  glColor3ub(255,255,0);
}

void DrawTextStats(GLWindow *glw, float tx, float ty)
{
  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // Draw
  w(); glRasterPos2f(  0.0f,  20.0f); printGLf(glw->font,"Speed:");
  y(); glRasterPos2f(120.0f,  20.0f); printGLf(glw->font,"%d",Speed);

  w(); glRasterPos2f(  0.0f,  40.0f); printGLf(glw->font,"Years:");
  y(); glRasterPos2f(120.0f,  40.0f); printGLf(glw->font,"%.2lf (%d/%llu)",((double)(Gen[Frame-1].g)-14)/((double)(MATINGCYCLES)), Frame-1, nGen-1);

  w(); glRasterPos2f(  0.0f,  60.0f); printGLf(glw->font,"Mating cycles/year:");
  y(); glRasterPos2f(120.0f,  60.0f); printGLf(glw->font,"%llu",MATINGCYCLES);
  
#if HOMINGDRIVE || HOMER || HOMER_KO || HOMER_HS
  w(); glRasterPos2f(  0.0f,  80.0f); printGLf(glw->font,"wildtype:");
  y(); glRasterPos2f(120.0f,  80.0f); printGLf(glw->font,"%llu",(Gen[Frame-1].geno_ze + Gen[Frame-1].geno_fo));
    
  w(); glRasterPos2f(  0.0f, 100.0f); printGLf(glw->font,"drive:");
  y(); glRasterPos2f(120.0f, 100.0f); printGLf(glw->font,"%llu",(Gen[Frame-1].geno_on + Gen[Frame-1].geno_fi));

#else
  w(); glRasterPos2f(  0.0f,  80.0f); printGLf(glw->font,"drive:");
  y(); glRasterPos2f(120.0f,  80.0f); printGLf(glw->font,"%llu",Gen[Frame-1].geno_tw);
      
  w(); glRasterPos2f(  0.0f, 100.0f); printGLf(glw->font,"wildtype:");
  y(); glRasterPos2f(120.0f, 100.0f); printGLf(glw->font,"%llu",(Gen[Frame-1].geno_on));

#endif
    
  w(); glRasterPos2f(  0.0f, 120.0f); printGLf(glw->font,"N/T/P(inoc):");
  y(); glRasterPos2f(120.0f, 120.0f); printGLf(glw->font,"%llu/%llu/%llu",INOC_INDIVIDUALS, NoInoc, InocPatchNo);
        
  w(); glRasterPos2f(  0.0f, 140.0f); printGLf(glw->font,"N(max):");
  y(); glRasterPos2f(120.0f, 140.0f); printGLf(glw->font,"%llu",GlobalMax);
    
  w(); glRasterPos2f(  0.0f, 160.0f); printGLf(glw->font,"T(reduction):");
  y(); glRasterPos2f(120.0f, 160.0f); printGLf(glw->font,"%.2lf",((double)(TPopSizeReduction-14))/((double)(MATINGCYCLES)));

  w(); glRasterPos2f(  0.0f, 180.0f); printGLf(glw->font,"Multiple Paternity:");
  y(); glRasterPos2f(120.0f, 180.0f); printGLf(glw->font,"%.2lf",((double)(Gen[Frame-1].mpat))/((double)(Gen[Frame-1].spat+Gen[Frame-1].mpat)));

}


void DrawSpeciesLegend(GLWindow *glw, float tx, float ty)
{
  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // Draw Column 1
  w(); glRasterPos2f(  0.0f,  20.0f); printGLf(glw->font,"Specialists:");
  glColor3ub(255,   0,   0);  glRasterPos2f(  0.0f,  40.0f); printGLf(glw->font,"Species 1");
  glColor3ub(255,   0, 255);  glRasterPos2f(  0.0f,  60.0f); printGLf(glw->font,"Species 2");
  glColor3ub(100,   0, 100);  glRasterPos2f(  0.0f,  80.0f); printGLf(glw->font,"Species 3");
  glColor3ub(255, 100,   0);  glRasterPos2f(  0.0f, 100.0f); printGLf(glw->font,"Species 4");

  w(); glRasterPos2f(  0.0f, 130.0f); printGLf(glw->font,"Generalists (2-r):");
  glColor3ub(  0,   0, 255);  glRasterPos2f(  0.0f, 150.0f); printGLf(glw->font,"Species 5");
  glColor3ub(  0, 100, 255);  glRasterPos2f(  0.0f, 170.0f); printGLf(glw->font,"Species 6");
  glColor3ub( 75,  75, 255);  glRasterPos2f(  0.0f, 190.0f); printGLf(glw->font,"Species 7");
  glColor3ub(  0, 180, 255);  glRasterPos2f(  0.0f, 210.0f); printGLf(glw->font,"Species 8");
  glColor3ub(  0, 255, 255);  glRasterPos2f(  0.0f, 230.0f); printGLf(glw->font,"Species 9");
  glColor3ub(100,   0, 255);  glRasterPos2f(  0.0f, 250.0f); printGLf(glw->font,"Species 10");

  // Draw Column 2
  w(); glRasterPos2f(120.0f,  20.0f); printGLf(glw->font,"Generalists (3-r):");
  glColor3ub(  0, 100,   0);  glRasterPos2f(120.0f,  40.0f); printGLf(glw->font,"Species 11");
  glColor3ub(  0, 200,  55);  glRasterPos2f(120.0f,  60.0f); printGLf(glw->font,"Species 12");
  glColor3ub(  0, 255,   0);  glRasterPos2f(120.0f,  80.0f); printGLf(glw->font,"Species 13");
  glColor3ub( 50, 130,   0);  glRasterPos2f(120.0f, 100.0f); printGLf(glw->font,"Species 14");

  w(); glRasterPos2f(120.0f, 130.0f); printGLf(glw->font,"Generalists (4-r):");
  glColor3ub(255, 255, 255);  glRasterPos2f(120.0f, 150.0f); printGLf(glw->font,"Species 15");

}


void DrawDemeGraph(GLWindow *glw, float tx, float ty)
{
  int i,j;  double t,c;

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // Draw the background patch character
  glBegin(GL_QUADS);
  for(i=0; i<Width; i++) {
    for(j=0; j<Height; j++) {

      t = ((double)Gen[Frame-1].p[i][j].n)/(NICHES-1);
#if BLACK_WHITE
      if( ((int)((double)(t)*((double)(NICHES-1))) == 0) ) glColor3ub(232,220,232);
      if( ((int)((double)(t)*((double)(NICHES-1))) == 1) ) glColor3ub(202,190,202);
      if( ((int)((double)(t)*((double)(NICHES-1))) == 2) ) glColor3ub(172,160,172);
      if( ((int)((double)(t)*((double)(NICHES-1))) == 3) ) glColor3ub(142,130,142);
        
#else
      glColor3ub(toByte(red  (t)*.3), toByte(green(t)*.3), toByte(blue (t)*.3));
#endif
      // t = ((double)Gen[Frame-1].p[i][j].n);//this is old
      // printf("\n %f", t);
#if PRESENTATION_MODE
      glColor3ub(toByte(red  (t)), toByte(green(t)), toByte(blue (t)));
#endif
      glVertex2i(i*Scale,           j*Scale);
      glVertex2i(i*Scale,           (j+1)*Scale);
      glVertex2i(i*Scale+Scale,     (j+1)*Scale);
      glVertex2i(i*Scale+Scale,     j*Scale);

#if PRETTY_PS
      if(PSFile1) {
        PS_Color(PSFile1, toByte(red(t)*.3)/255., toByte(green(t)*.3)/255., toByte(blue(t)*.3)/255.);
        PS_Box(PSFile1, i/2., j/2., .5, .5, 1);
      }
#endif

    }
  }
  glEnd();

  // Draw the patch size circle and outline
  for(i=0; i<Width && Gen[Frame-1].Popsize; i++) {
    for(j=0; j<Height; j++) {
      // Skip empty patches
      if( !Gen[Frame-1].p[i][j].ni ) continue;

#if SCALING
      t = Gen[Frame-1].p[i][j].ni/((double)Gen[Frame-1].rMax);
#else
      //t = Gen[Frame-1].p[i][j].ni/((Gen[Frame-1].rMax<Max/_b)?(Max/_b):((double)Gen[Frame-1].rMax));
      t = ((double)(Gen[Frame-1].p[i][j].ni))/((double)(PPopMax));
      
#endif

#if BLACK_WHITE
      c = Gen[Frame-1].p[i][j].sid;
      if( t >= 0 ){

      //  "#cd34b5", female semi 1
      //  "#9d02d7", female drive 2
      //  "#fa8775", male drive 3
      //  "#FFE37E", male semi 4
#if HOMER_HS
      if( (((int)(c)) == 0) )  glColor3ub(255, 255, 255);
      if( (((int)(c)) == 1) )  glColor3ub(240, 240, 240);
      if( (((int)(c)) == 2) )  glColor3ub(255,   0,   0);
      if( (((int)(c)) == 3) )  glColor3ub(152,   0,  51);
      if( (((int)(c)) == 4) )  glColor3ub(255, 250,   0);
      if( (((int)(c)) == 5) )  glColor3ub(255, 200,   0);
      if( (((int)(c)) == 6) )  glColor3ub(183,   0, 255);
      if( (((int)(c)) == 7) )  glColor3ub(123,   0, 255);
#else
          
      if( (((int)(c)) == 0) )  glColor3ub( 80,  80,  80);
      if( (((int)(c)) == 5) )  glColor3ub( 80,  80,  80);
      if( (((int)(c)) == 1) )  glColor3ub(255,   0,   0);
      if( (((int)(c)) == 2) )  glColor3ub(255,   0,   0);
      if( (((int)(c)) == 3) )  glColor3ub(153,  51, 255);
      if( (((int)(c)) == 4) )  glColor3ub(  0, 153, 255);
      if( (((int)(c)) == 5) )  glColor3ub(  0, 153, 255);
      if( (((int)(c)) == 6) )  glColor3ub(102,   0, 204);
      if( (((int)(c)) == 7) )  glColor3ub(  0,  51, 204);
      if( (((int)(c)) == 8) )  glColor3ub(255, 100,   0);
      if( (((int)(c)) == 9) )  glColor3ub(  0, 180, 255);
      if( (((int)(c)) == 10) ) glColor3ub(  0, 255, 255);
      if( (((int)(c)) == 11) ) glColor3ub(  0, 200,  55);
      if( (((int)(c)) == 12) ) glColor3ub(100,   0, 255);
      if( (((int)(c)) == 13) ) glColor3ub(  0, 255,   0);
      if( (((int)(c)) == 14) ) glColor3ub( 50, 130,   0);
      if( (((int)(c)) == 15) ) glColor3ub(255, 255, 255);
      if( (((int)(c)) > 15) )  glColor3ub( 50,  50,  50);
#endif
      glCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,sqrt((Scale*.5)*(Scale*.5)*t));
      glOutlineCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,sqrt((Scale*.5)*(Scale*.5)*t));
      }

#else
      if( ColorMode == COLOR_PREF ) {
        // Circle colors are based on niche preference
        c = Gen[Frame-1].p[i][j].pnf/((double)(NICHES-1));
#if PRESENTATION_MODE
        glColor3ub(0,0,0);
        glCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,(Scale*sqrt(t))*.5);
        glColor3ub(toByte(red  (c)), toByte(green(c)), toByte(blue (c)));      
        if((Scale*t)*.5-2 > 1)
          glCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,(Scale*sqrt(t))*.5-2);
#else
        glColor3ub(toByte(red  (c)*.3), toByte(green(c)*.3), toByte(blue (c)*.3));
        if( t >= 0 )
          glCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,sqrt((Scale*.5)*(Scale*.5)*t));
        glColor3ub(toByte(red  (c)), toByte(green(c)), toByte(blue (c)));
        if( t >= 0 )
          glOutlineCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,sqrt((Scale*.5)*(Scale*.5)*t));      
#endif

#if PRETTY_PS
        if(PSFile1) {
          PS_Color(PSFile1, toByte(red(c)*.3)/255., toByte(green(c)*.3)/255., toByte(blue(c)*.3)/255.);
          if( t >= 0 )
            PS_Circle(PSFile1, i/2.+.25, j/2.+.25, sqrt(t)*.25, 1);
          PS_Color(PSFile1, toByte(red(c))/255., toByte(green(c))/255., toByte(blue(c))/255.);
          if( t >= 0 )
            PS_Circle(PSFile1, i/2.+.25, j/2.+.25, sqrt(t)*.25, 0);
        }
#endif
      } else {
        // Circle color is based on marker trait, outline is still preference
        if( !Filter[Gen[Frame-1].p[i][j].n] ) {
          // Only draw if not filtered
          c = Gen[Frame-1].p[i][j].k;
          glColor3f(c, c, c);
          if( t >= 0 )
            glCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,sqrt((Scale*.5)*(Scale*.5)*t));
          c = Gen[Frame-1].p[i][j].pnf/((double)(NICHES-1));
          glColor3ub(toByte(red  (c)), toByte(green(c)), toByte(blue (c)));
          if( t >= 0 )
            glOutlineCircle(i*Scale+Scale*.5,j*Scale+Scale*.5,sqrt((Scale*.5)*(Scale*.5)*t));
        }
      }
#endif
    }
  }

  // Outline
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(0,0);
  glVertex2i(Width*Scale,0);
  glVertex2i(Width*Scale,Height*Scale);
  glVertex2i(0,Height*Scale);
  glEnd();
}


////////////////////////////////////////////////////////////////////////////////
// Within-Deme graph
//
// This gives a surface graph of a property of a single deme, where the x-axis
// and the z-axis are environment height and width, respectively.
////////////////////////////////////////////////////////////////////////////////

void Draw3DDeme(GLWindow *glw)
{
  u32b x,y;
  //char title[64];
  patch *p;
  int i;

  ////////////////////////////////////////////////////////////
  // 3D trait cluster view
  ////////////////////////////////////////////////////////////

  // Outline widget border
  Yellow();
  glBegin(GL_LINE_LOOP);
  glVertex2f(0.0f,0.0f);
  glVertex2f(0.0f,1.0f);
  glVertex2f(1.0f,1.0f);
  glVertex2f(1.0f,0.0f);
  glEnd();

  // Draw title
  //glRasterPos2f(0.5f-((strlen("traits")*6.0f/2.0f)/Plot.w),0.9f);
  //printGLf(glw->font,"traits");
  //sprintf(title,"Min: %d",Plot.spctrshld);
  //glRasterPos2f(0.5f-((strlen(title)*6.0f/2.0f)/Plot.w),0.95f);
  //printGLf(glw->font,title);

  // Switch over to a custom viewport "mode" which will let
  // us draw with our own perspective within our widget boundries.
  ViewPort3D(Plot.x, glw->height-Plot.y-Plot.h, Plot.w, Plot.h);

  // Rotate and translate to position 0.5 at the origin.
  // This allows all further drawing to be done [0,1], which,
  // IMHO, is very nice to work with.
  // glTranslatef(x,yf,z);                                        // Move graph "object" up/down on screen
  glRotatef(-55.0f+Plot.rot.x+(Plot.ny-Plot.oy),1.0f,0.0f,0.0f);  // Tilt back (rotate on graph's x axis)
  glRotatef(Plot.rot.z+(Plot.nx-Plot.ox),0.0f,0.0f,1.0f);         // Rotate on graph's y-axis
  glTranslatef(-0.5f,-0.5f,-0.0f);                                // Position .5 at origin
  glScalef(1.0f,1.0f,0.75f);

  // Outline bottom and top of heightmap
  glColor3f(0.3f, 0.3f, 0.3f);
  glEnable(GL_LINE_SMOOTH);
  glBegin(GL_LINE_LOOP); 
  glVertex2f(0.0f,0.0f);
  glVertex2f(0.0f,1.0f);
  glVertex2f(1.0f,1.0f);
  glVertex2f(1.0f,0.0f);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f(0.0f,0.0f,1.0f);
  glVertex3f(0.0f,1.0f,1.0f);
  glVertex3f(1.0f,1.0f,1.0f);
  glVertex3f(1.0f,0.0f,1.0f);
  glEnd();
  glDisable(GL_LINE_SMOOTH);

  //
  // 3-d plot
  //
  // Draw the plot
  glColor3f(1.0f, 0.0f, 1.0f);
  glPointSize(3.0f);
  glEnable(GL_POINT_SMOOTH);
  glBegin(GL_POINTS);
  for(i=0; i<nGen; i++) {
    for(x=0; x<Width; x++) {
      for(y=0; y<Height; y++) {
        p = &Gen[Frame-1].p[x][y];
        // Color based on fourth trait
        glColor3f(0.0f, p->af[3], 1.0f-p->af[3]);
        // Position is based on the three other traits
        glVertex3f(p->af[0], p->af[1], p->af[2]);
      }
    }
  }

  glEnd();
  glPointSize(1.0f);
  glDisable(GL_POINT_SMOOTH);

  // Label the four corners
  White();
  glRasterPos3f(0.0f,0.0f,0.1f); printGLf(glw->font,"(0,0,0)");
  glRasterPos3f(0.0f,1.0f,0.1f); printGLf(glw->font,"(0,1,0)");
  glRasterPos3f(1.0f,1.0f,0.1f); printGLf(glw->font,"(1,1,0)");
  glRasterPos3f(1.0f,0.0f,0.1f); printGLf(glw->font,"(1,0,0)");

  // Restore 2d viewport for other widgets
  ViewPort2D(glw);
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#define HS (((Scale-HackScale)<0)?(0):(Scale-HackScale))

void DrawNicheHistogram(GLWindow *glw, float tx, float ty)
{
  int i; double t,c;

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // The outline for the histogram should be drawn first
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(-1,0);
  glVertex2i(NICHES*(20+HS)+1,0);
  glVertex2i(NICHES*(20+HS)+1,100);
  glVertex2i(-1,100);
  glEnd();

  // Draw our "niche histogram"
  for(i=0; i<NICHES && Gen[Frame-1].Popsize; i++) {
    t = ((double)i)/(NICHES-1);
#if BLACK_WHITE
    if( ((int)((double)(t)*((double)(NICHES-1))) == 0) ) glColor3ub(232,220,232);
    if( ((int)((double)(t)*((double)(NICHES-1))) == 1) ) glColor3ub(202,190,202);
    if( ((int)((double)(t)*((double)(NICHES-1))) == 2) ) glColor3ub(172,160,172);
    if( ((int)((double)(t)*((double)(NICHES-1))) == 3) ) glColor3ub(142,130,142);
#else
    glColor3ub(toByte(red  (t)*.3), toByte(green(t)*.3), toByte(blue (t)*.3));
#endif
#if PRESENTATION_MODE
    glColor3ub(toByte(red  (t)*.7), toByte(green(t)*.7), toByte(blue (t)*.7));

#endif
    glBegin(GL_QUADS);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].Nichef[i])/Gen[Frame-1].mNichef)*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].Nichef[i])/Gen[Frame-1].mNichef)*99);
    glEnd();    
    // Draw a tick mark to represent the proportion of "real" species
#if BLACK_WHITE
    glColor3ub(60,60,60);
#else
    glColor3ub(toByte(red  (t)), toByte(green(t)), toByte(blue (t)));
#endif

#if PRESENTATION_MODE
    glColor3ub(toByte(red  (1-t)), toByte(green(1-t)), toByte(blue (1-t)));

#endif
    c = (((double)Gen[Frame-1].Prefs[i][i]) / Gen[Frame-1].mPref) * 100;
    glBegin(GL_LINES);
    glVertex2i(i*(20+HS),100-c);
    glVertex2i((i+1)*(20+HS),100-c);
    glEnd();
    if( i == Gen[Frame-1].Mpnf ) { 
      // Draw a marker on the bar which represents the most prefered niche (fitness)
#if PRESENTATION_MODE
      glRasterPos2f(i*(20+HS)+20, 90.0f);
#else
      glRasterPos2f(i*(20+HS)+13, 90.0f);
#endif
      printGLf(glw->font,"fi");
    }
    if( i == Gen[Frame-1].Mpnm ) {
      // Draw a marker on the bar which represents the most prefered niche (migration)
#if PRESENTATION_MODE
      glRasterPos2f(i*(20+HS)+3, 90.0f);
#else
      glRasterPos2f(i*(20+HS)+3, 90.0f);
#endif
      printGLf(glw->font,"mi");
    }
    if( ColorMode != COLOR_PREF ) {
      if( !Filter[i] ) {
        glColor3ub(255,255,255);
        glCircle(i*(20+HS)+((20.0+HS)/2.),50,2);
      }
    }
    glLogicOp(GL_COPY);
  }
}

void DrawSpeciesHistogramI(GLWindow *glw, float tx, float ty)
{
  int i,sp=0; double t;
  sp = ((int)(pow((double)(2),(double)(nTraits))));

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // The outline for the histogram should be drawn first
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(-1,0);
  glVertex2i(sp*(20+HS)+1,0);
  glVertex2i(sp*(20+HS)+1,100);
  glVertex2i(-1,100);
  glEnd();

  // Draw our "species histogram"
  for(i=0; i<sp; i++) {
    t = ((double)(i))/((double)(sp));
    if( ((int)((double)(t)*((double)(sp))) == 0) )      glColor3ub( 50,  50,  50);
    if( ((int)((double)(t)*((double)(sp))) == 1) )      glColor3ub(255,   0,   0);
    if( ((int)((double)(t)*((double)(sp))) == 2) )      glColor3ub(255,   0, 255);
    if( ((int)((double)(t)*((double)(sp))) == 3) )      glColor3ub(  0,   0, 255);
    if( ((int)((double)(t)*((double)(sp))) == 4) )      glColor3ub(100,   0, 100);
    if( ((int)((double)(t)*((double)(sp))) == 5) )      glColor3ub(  0, 100, 255);
    if( ((int)((double)(t)*((double)(sp))) == 6) )      glColor3ub( 75,  75, 255);
    if( ((int)((double)(t)*((double)(sp))) == 7) )      glColor3ub(  0, 100,   0);
    if( ((int)((double)(t)*((double)(sp))) == 8) )      glColor3ub(255, 100,   0);
    if( ((int)((double)(t)*((double)(sp))) == 9) )      glColor3ub(  0, 180, 255);
    if( ((int)((double)(t)*((double)(sp))) == 10) )     glColor3ub(  0, 255, 255);
    if( ((int)((double)(t)*((double)(sp))) == 11) )     glColor3ub(  0, 200,  55);
    if( ((int)((double)(t)*((double)(sp))) == 12) )     glColor3ub(100,   0, 255);
    if( ((int)((double)(t)*((double)(sp))) == 13) )     glColor3ub(  0, 255,   0);
    if( ((int)((double)(t)*((double)(sp))) == 14) )     glColor3ub( 50, 130,   0);
    if( ((int)((double)(t)*((double)(sp))) == 15) )     glColor3ub(255, 255, 255);

    if( ((double)Gen[Frame-1].sphist.sphsize[i]) == 0. || ((double)Gen[Frame-1].sphist.msphsize) == 0. ){
      glBegin(GL_QUADS);
      glVertex2i(i*(20+HS),100);
      glVertex2i(i*(20+HS),99);
      glVertex2i((i+1)*(20+HS),99);
      glVertex2i((i+1)*(20+HS),100);
      glEnd();
    } else {
      glBegin(GL_QUADS);
      glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].sphist.sphsize[i])/((double)Gen[Frame-1].sphist.msphsize))*99);
      glVertex2i(i*(20+HS),99);
      glVertex2i((i+1)*(20+HS),99);
      glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].sphist.sphsize[i])/((double)Gen[Frame-1].sphist.msphsize))*99);
      glEnd();
    }
    glLogicOp(GL_COPY);
  }
}

/*
void DrawSpeciesHistogramII(GLWindow *glw, float tx, float ty)
{
  int sp, spo, s;
  sp=spo=s=0;
  sp = ((int)(pow((double)(2),(double)(nTraits))));
  spo= (2*((int)(pow((double)(2),(double)(nTraits)))));
  s  = spo - sp;


  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // The outline for the histogram should be drawn first
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(-1,0);
  glVertex2i(sp*(20+HS)+1,0);
  glVertex2i(sp*(20+HS)+1,100);
  glVertex2i(-1,100);
  glEnd();

  // Draw our "species histogram"
  int j;
  for(j=0; j<sp; j++) {
    glColor3ub(100, 100, 100);
    glBegin(GL_QUADS);
    glVertex2i((j)*(20+HS),100-(((double)Gen[Frame-1].sphist.sphsize[j+s])/((double)Gen[Frame-1].sphist.msphsize))*99);
    glVertex2i((j)*(20+HS),99);
    glVertex2i((j+1)*(20+HS),99);
    glVertex2i((j+1)*(20+HS),100-(((double)Gen[Frame-1].sphist.sphsize[j+s])/((double)Gen[Frame-1].sphist.msphsize))*99);
    glEnd();

    glLogicOp(GL_COPY);
  }
}
*/

void DrawPreferenceGraph(GLWindow *glw, float tx, float ty)
{
  int i,j;  double t;

#if PRETTY_PS
  char b[128];
#endif

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

#if PRETTY_PS
  if(PSFile0 && Gen[Frame-1].Popsize) {
    // Move the graph over so there is room for labels
    PS_Translate(PSFile0,.25,0);
    // Draw a white background 
    PS_Gray(PSFile0, 1.0);
    PS_Box(PSFile0,0,0,NICHES*.5,NICHES*.5,1);
  }
#endif

  /* Draw preference graph */
  glBegin(GL_QUADS);
  for(i=0; i<NICHES && Gen[Frame-1].Popsize; i++)
    for(j=0; j<NICHES; j++) {
      t = (((double)Gen[Frame-1].Prefs[i][j]) / Gen[Frame-1].mPref) * 255;
      glColor3ub(toByte(t), toByte(t), toByte(t));
      glVertex2i(i*(20+HS),j*(20+HS));
      glVertex2i(i*(20+HS),(j+1)*(20+HS));
      glVertex2i((i+1)*(20+HS),(j+1)*(20+HS));
      glVertex2i((i+1)*(20+HS),j*(20+HS));
#if PRETTY_PS
      if(PSFile0) {
        // Draw the boxes
        t = toByte(t)/255.;
#if LOG_SCALING
        t = log(toByte(t)+1)/log(256);
#endif
        t = 1.-t;
        PS_Gray(PSFile0,t);
        PS_Box(PSFile0,i/2.,j/2.,.5,.5,1);
      }
#endif
    }
  glEnd();

#if PRETTY_PS
  if(PSFile0 && Gen[Frame-1].Popsize) {
    // Draw an outline
    PS_Gray(PSFile0,0.);
    PS_Box(PSFile0,0,0,NICHES*.5,NICHES*.5,0); 
    // Restore the origin and draw some text labels
    PS_Translate(PSFile0,-.25,0);
    // Draw left labels
    for(i=0; i<NICHES; i++) {
      sprintf(b,"%d",i);
      PS_Text(PSFile0, .115, i*.5+.25, 12, b);
    }
    // Draw bottom labels
    for(i=0; i<NICHES; i++) {
      sprintf(b,"%d",i);
      PS_Text(PSFile0, i*.5+.45, NICHES*.5+.16, 12, b);
    }
  }
#endif

  /* Outline */
  glBegin(GL_LINES);
  for(i=0; i<NICHES;) {
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i(i*(20+HS),0);
    i++;
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i(i*(20+HS),0);
    i--;
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i(i*(20+HS),(20+HS)*NICHES);
    i++;
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i(i*(20+HS),(20+HS)*NICHES);
    i--;
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i(0,i*(20+HS));
    i++;
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i(0.0,i*(20+HS));
    i--;
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i((20+HS)*NICHES,i*(20+HS));
    i++;
    t = ((double)i)/(NICHES-1);
    glColor3ub(toByte(red(t)), toByte(green(t)), toByte(blue(t)));
    glVertex2i((20+HS)*NICHES,i*(20+HS));
  }
  glEnd();
}


void ProcessDuration(int spc, int start, int end, int flag)
{
  static double  mdur,adur,amrange,count,amrangemaxdur;
  static FILE   *f, *g;
  static int     mdurs, mdure;
  int            i;

  switch(flag) {
    case -1:
      // per-file
      if( !(f=fopen("durationdata","w")) ) {
        fprintf(stderr, "Error: cannot open duration data file!\n");
      }
      if( !(g=fopen("maxdurationdata","w")) ) {
        fprintf(stderr, "Error: cannot open max duration data file!\n");
      }
      break;

    case 0:
      // init (per species)
      count           = 0.0; // number of separate durations
      adur            = 0.0; // average duration
      mdur            = 0.0; // maximum duration
      mdurs           = 0;   // maximum duration start
      mdure           = 0;   // maximum duration end
      amrange         = 0.0; // summing all the ranges across duration specified
      amrangemaxdur   = 0.0; // average maximum range during the maximum duration
      break;

    case 1:
      // Still causing "duration" to be applied to average (per species)
      // printf("duration: spc%d [%d,%d) %d\n", Gen[start].range.idspecies[spc], start, end, end-start);
      adur += end-start;                          // summing the durations to calculate average length of durations
      count++;                                    // number of durations
      if( end-start > mdur ) {                    // find max duration length meanwhile
        mdur = end-start;
        mdurs = start;                            // set the start of max duration
        mdure = end;                              // set the end of max duration
      }

      // Find average max range for the current duration
      for(i=start; i<end; i++) {
        amrange += Gen[i].range.mrange[spc];      //summing all ranges for average range calculation later
      }
      break;

    case 2:
      // done with all "durations" collected so far (per species)
      if( count ){
        amrange /= adur;  // doesn't seem very useful
        adur    /= count; // doesn't seem very useful
        // also calculate the mean max range during the max duration
        for(i=mdurs; i<mdure; i++){
          amrangemaxdur += Gen[i].range.mrange[spc];      //summing all ranges for average range calculation later
        }
        amrangemaxdur /= mdur;
      } else {
        adur =  0.;
      }

      // printf("duration: spc%d has avg duration len (with %lf samples) of %lf\n", spc, count, adur);
      // for "durationdata" file
      fprintf(f,"%lf ", adur);
      for(i=0; i<spc; i++) { 
        fprintf(f,"0.000000 ");
      }
      fprintf(f,"%lf ",amrange);
      for(i=spc+1; i<(((int)(pow((double)(2),(double)(nTraits))))); i++) { 
        fprintf(f,"0.000000 ");
      }
      fprintf(f,"\n");

      // for "maxdurationdata" file
      fprintf(g,"%lf ", mdur);
      for(i=0; i<spc; i++) {
        fprintf(g,"0.000000 ");
      }
      fprintf(g,"%lf ",amrangemaxdur);

      printf("species%dmaxdur: %lff-species%drangemdur: %lff\n", spc, mdur, spc, amrangemaxdur);

      for(i=spc+1; i<(((int)(pow((double)(2),(double)(nTraits))))); i++) {
        fprintf(g,"0.000000 ");
      }
      fprintf(g,"\n");
      break;

    case 3:
      // per file
      fclose(f);
      fclose(g);
      break;
  }
}


void RecordInvasion()
{
  int     i, l;

  for( i=0; i<nGen; i++){
    // time to colonization
    if( Timetocolonize == 0 ){
      for(l=0; l<Width; l++){
        if( Gen[i].p[63][l].ni  > 0 ){
          Timetocolonize = i;
          break;
        }
      }
    }
  }
}

void RecordReduction()
{
  int     i;
  double  r;
  
  for( i = FirstInocGen; i<nGen; i++){
    // time to colonization
    r = (((double)(Gen[i].Popsize))/((double)(GlobalMax)));
    // printf("r: %lf, i: %llu\n", r, i);
    if( r < 0.01 ){
      TPopSizeReduction = i;
      //printf("%llu\n", TPopSizeReduction);
      break;
    }
  }
}


void WriteDemes()
{
  int i,j;  double t,c;
  FILE *f;  double sx=4.0/(Width*Scale),sy=4.0/(Height*Scale);

  // Open eps file
  if( !(f = PS_Start("demesfig.eps",Width*Scale*sx,Height*Scale*sy)) ) {
    fprintf(stderr,"Could not open output post-script file! (demesfig.eps)\n");
    exit(1);
  }
//printf( "demesfig.eps file is opened!\n" );
  // Move the graph over so there is room for labels
  PS_Translate(f,0.,0.);

  // Draw a white background
  PS_Color(f, 1.0, 1.0, 1.0);
  PS_Quad(f,
          0,              0,
          0,              sy*Scale*Height,
          sx*Scale*Width, sy*Scale*Height,
          sx*Scale*Width, 0,
          1);


  // Draw the background patch character
  for(i=0; i<Width; i++) {
    for(j=0; j<Height; j++) {

      t = ((double)Gen[Frame-1].p[i][j].n)/(NICHES-1);
#if BLACK_WHITE
      if( ((int)((double)(t)*((double)(NICHES-1))) == 0) ) PS_Color(f,232/((float)255),220/((float)255),232/((float)255));
      if( ((int)((double)(t)*((double)(NICHES-1))) == 1) ) PS_Color(f,202/((float)255),190/((float)255),202/((float)255));
      if( ((int)((double)(t)*((double)(NICHES-1))) == 2) ) PS_Color(f,172/((float)255),160/((float)255),172/((float)255));
      if( ((int)((double)(t)*((double)(NICHES-1))) == 3) ) PS_Color(f,142/((float)255),130/((float)255),142/((float)255));
#else
      PS_Color(f,toByte(red  (t)*.3)/((float)255), toByte(green(t)*.3)/((float)255), toByte(blue (t)*.3)/((float)255));
#endif
      //t = ((double)Gen[Frame-1].p[i][j].n);//this is old
      //printf("\n %f", t);
      PS_Quad(f,
              i*Scale*sx,           j*Scale*sy,
              i*Scale*sx,           (j+1)*Scale*sy,
              (i*Scale+Scale)*sx,   (j+1)*Scale*sy,
              (i*Scale+Scale)*sx,   j*Scale*sy,
              1);
    }
  }

  // Draw the patch size circle and outline 
  for(i=0; i<Width && Gen[Frame-1].Popsize; i++) {
    for(j=0; j<Height; j++) {
      /* Skip empty patcehs */
      if( !Gen[Frame-1].p[i][j].ni )
        continue;
      //t = Gen[Frame-1].p[i][j].ni/((Gen[Frame-1].rMax<Max/_b)?(Max/_b):((double)Gen[Frame-1].rMax));

        t = Gen[Frame-1].p[i][j].ni/((double)(PPopMax));
        
        
#if BLACK_WHITE
      c = Gen[Frame-1].p[i][j].sid;
      if( t >= 0 ){
        if( (((int)(c)) == 0) )  PS_Color(f, 100/((float)255), 100/((float)255), 100/((float)255));
        if( (((int)(c)) == 1) )  PS_Color(f, 255/((float)255),   0/((float)255),   0/((float)255));
        if( (((int)(c)) == 2) )  PS_Color(f, 186/((float)255),  85/((float)255), 211/((float)255));
        if( (((int)(c)) == 3) )  PS_Color(f, 153/((float)255),  51/((float)255), 255/((float)255));
        if( (((int)(c)) == 4) )  PS_Color(f,  80/((float)255),  80/((float)255),  80/((float)255));
        if( (((int)(c)) == 5) )  PS_Color(f, 152/((float)255),   0/((float)255),  51/((float)255));
        if( (((int)(c)) == 6) )  PS_Color(f, 102/((float)255),   0/((float)255), 204/((float)255));
        if( (((int)(c)) == 7) )  PS_Color(f,   0/((float)255),  51/((float)255), 204/((float)255));
        if( (((int)(c)) == 8) )  PS_Color(f, 255/((float)255), 100/((float)255),   0/((float)255));
        if( (((int)(c)) == 9) )  PS_Color(f,   0/((float)255), 180/((float)255), 255/((float)255));
        if( (((int)(c)) == 10) ) PS_Color(f,   0/((float)255), 255/((float)255), 255/((float)255));
        if( (((int)(c)) == 11) ) PS_Color(f,   0/((float)255), 200/((float)255),  55/((float)255));
        if( (((int)(c)) == 12) ) PS_Color(f, 100/((float)255),   0/((float)255), 255/((float)255));
        if( (((int)(c)) == 13) ) PS_Color(f,   0/((float)255), 255/((float)255),   0/((float)255));
        if( (((int)(c)) == 14) ) PS_Color(f,  50/((float)255), 130/((float)255),   0/((float)255));
        if( (((int)(c)) == 15) ) PS_Color(f, 255/((float)255), 255/((float)255), 255/((float)255));
        if( (((int)(c)) >  15) ) PS_Color(f,  50/((float)255),  50/((float)255),  50/((float)255));
        PS_Circle(f, (i*Scale+Scale*.5)*sx, (j*Scale+Scale*.5)*sy, sqrt((Scale*.5)*(Scale*.5)*t)*sy, 1);
        PS_Circle(f, (i*Scale+Scale*.5)*sx, (j*Scale+Scale*.5)*sy, sqrt((Scale*.5)*(Scale*.5)*t)*sy, 0);
      }
#else
      // Circle colors are based on niche preference
      c = Gen[Frame-1].p[i][j].pnf/((double)(NICHES-1));
      PS_Color(f,toByte(red  (c)*.3)/((float)255), toByte(green(c)*.3)/((float)255), toByte(blue (c)*.3)/((float)255));
      if( t >= 0 ) {
        PS_Circle(f, (i*Scale+Scale*.5)*sx, (j*Scale+Scale*.5)*sy, sqrt((Scale*.5)*(Scale*.5)*t)*sy, 1);
      }
      PS_Color(f, toByte(red  (c))/((float)255), toByte(green(c))/((float)255), toByte(blue (c))/((float)255));
      if( t >= 0 ) {
        PS_Circle(f, (i*Scale+Scale*.5)*sx, (j*Scale+Scale*.5)*sy, sqrt((Scale*.5)*(Scale*.5)*t)*sy, 0);
      }
#endif
    }
  }

  // Outline
  PS_Color(f,0.f,0.f,0.f);
  PS_Box(f, 0.f, 0.f, Width*Scale*sx, Height*Scale*sy, 0);

  // Restore origin and write text labels
  PS_Translate(f,-0.5,0.);
  // Close eps file
  PS_End(f);
}


#if READ_DISTANCE
void DrawDVD(GLWindow *glw, float tx, float ty)
{
  int i,j;  double t,x,y;

//#if PRETTY_PS
  // Write to disk if eps file is open
  //if( PSFile0 ) {
  //  WriteDVD();
  //}
//#endif

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // Draw dvd
  x = (8*(20+HS))/((double)12);
  y = (8*(20+HS))/((double)DVDy);
  glBegin(GL_QUADS);
  for(i=0; i<12; i++) {
    for(j=0; j<DVDy; j++) {
      t = (((double)Gen[Frame-1].dvd[i][j]) / Gen[Frame-1].mdvd) * 255;
      glColor3ub(toByte(t), toByte(t), toByte(t));
      glVertex2i(  i  *x,   j  *y);
      glVertex2i(  i  *x, (j+1)*y);
      glVertex2i((i+1)*x, (j+1)*y);
      glVertex2i((i+1)*x,   j  *y);
    }
  }
  glEnd();

  // Outline and labels
  glBegin(GL_LINE_LOOP);
  glColor3ub(255,255,255);
  glVertex2i(0,            0);
  glVertex2i(0,            8*(20+HS));
  glVertex2i(8*(20+HS), 8*(20+HS));
  glVertex2i(8*(20+HS), 0);
  glEnd();
  glRasterPos2f(4*(5+HS)-16, 8*(20+HS)+12.0f);   printGLf(glw->font,"0        density      1");
  glRasterPos2f(-12,          4*(20+HS)-60-5);   printGLf(glw->font,"1");
  glRasterPos2f(-12,          4*(20+HS)-20-5);   printGLf(glw->font,"d");
  glRasterPos2f(-12,          4*(20+HS)-10-5);   printGLf(glw->font,"i");
  glRasterPos2f(-12,          4*(20+HS)-5);      printGLf(glw->font,"s");
  glRasterPos2f(-12,          4*(20+HS)+10-5);   printGLf(glw->font,"t");
  glRasterPos2f(-12,          4*(20+HS)+20-5);   printGLf(glw->font,"a");
  glRasterPos2f(-12,          4*(20+HS)+30-5);   printGLf(glw->font,"n");
  glRasterPos2f(-12,          4*(20+HS)+40-5);   printGLf(glw->font,"c");
  glRasterPos2f(-12,          4*(20+HS)+50-5);   printGLf(glw->font,"e");
  glRasterPos2f(-12,          4*(20+HS)+80-5);   printGLf(glw->font,"0");
}
#endif


void DrawGenerationSlider(GLWindow *glw, float tx, float ty)
{
  double t;

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  /* Draw generation slider */
  glColor3ub(255, 255, 255);
  glBegin(GL_LINES);
  /* Bar */
  glVertex2i(0, 15.0f);
  glVertex2i(Width*Scale, 15.0f);
  /* Vertical "tick" mark */
  t = (((double)(Frame-1))/(nGen-1))*(Width*Scale);
  glVertex2i(t, 0.0f);
  glVertex2i(t, 30.0f);
  glEnd();
}

void DrawSummaryGraph(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0; static double *buckets[6];

  if ( nGen <= 1 )
    return;

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
    // Draw graph here from buckets
    glBegin(GL_LINES);
    
    // gcc
    glColor3ub(255,   0,   0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
    }
    // amk
    glColor3ub(250, 250, 250);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
    }
    // amc
    glColor3ub(100, 100, 100);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
    }
    // 
    glColor3ub(255,   0,   0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[4][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[4][i]);
    }

    glEnd();
    // Outline
    glColor3ub(255, 255, 255);
    glBegin(GL_LINE_LOOP);
    glVertex2i(0,0);
    glVertex2i(Width*Scale,0);
    glVertex2i(Width*Scale,150);
    glVertex2i(0,150);
    glEnd();
    // Labels for the summary graph
    i = -40;

    glColor3ub(100, 100, 100);  glRasterPos2f(i+=40, 165);  printGLf(glw->font,"Wildtype (L1)");
    glColor3ub(255,  0,    0);  glRasterPos2f(i+=85, 165);  printGLf(glw->font,"Drive (L1)");
    glColor3ub(250, 250, 250);  glRasterPos2f(i+=95, 165);  printGLf(glw->font,"Infertile (L2)");
    
    
  } else {
    // Initialize buckets
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<6; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    // Fill in buckets
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));
      buckets[0][j] += Gen[i].acp;
      buckets[1][j] += Gen[i].amk;
      buckets[2][j] += Gen[i].amc;
      buckets[3][j] += 0.;
      buckets[4][j] += 0.;
      // Just count; used for scaling
      buckets[5][j]++;
    }
    // Scale buckets
    for(i=0; i<6; i++)
      for(j=0; j<nb; j++) {
        buckets[i][j] = 150-((buckets[i][j]/buckets[5][j])*150);
      }

    // Jump tp drawing section to cause the graph to actually be drawn
    goto draw;
  }
}


void DrawGenotypeHistogram(GLWindow *glw, float tx, float ty)
{
  int i,sp=0; double t;
  sp = 8;

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  // The outline for the histogram should be drawn first
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(-1,0);
  glVertex2i(sp*(20+HS)+1,0);
  glVertex2i(sp*(20+HS)+1,100);
  glVertex2i(-1,100);
  glEnd();
    
    
  i = 0;
    // "/((double)Gen[Frame-1].Popsize))*99"
    glBegin(GL_QUADS);
    glColor3ub(200, 200, 200);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_ze)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_ze)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
    glBegin(GL_QUADS);
    glColor3ub(255,   0,   0);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_on)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_on)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
    glBegin(GL_QUADS);
    glColor3ub(186,  85, 211);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_tw)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_tw)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
    glBegin(GL_QUADS);
    glColor3ub(  0, 150,   0);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_th)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_th)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
    glBegin(GL_QUADS);
    glColor3ub(150, 150, 150);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_fo)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_fo)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
    glBegin(GL_QUADS);
    glColor3ub(200,   0,   0);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_fi)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_fi)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
    glBegin(GL_QUADS);
    glColor3ub(138,  43, 226);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_si)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_si)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
    glBegin(GL_QUADS);
    glColor3ub(  0, 100,   0);
    glVertex2i(i*(20+HS),100-(((double)Gen[Frame-1].geno_se)/((double)Gen[Frame-1].Popsize))*99);
    glVertex2i(i*(20+HS),99);
    glVertex2i((i+1)*(20+HS),99);
    glVertex2i((i+1)*(20+HS),100-(((double)Gen[Frame-1].geno_se)/((double)Gen[Frame-1].Popsize))*99);
    glEnd();
  i++;
        
    glLogicOp(GL_COPY);
}


void DrawGenotypePop(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0; static double *buckets[9];

  if ( nGen <= 1 )
    return;

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
    // Draw graph here from buckets
    glBegin(GL_LINES);

    //    0. wildtype female
    glColor3ub(200, 200, 200);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
    }
    
    //    1. drive female
    glColor3ub(255,   0,   0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
    }
    
    //    2. infertile female/or resistant (vs other in HOMER_HS)
    glColor3ub(186,  85, 211);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
    }

    //    3. other (vs rescue lost in HOMER_HS) female
    glColor3ub(  0, 150,   0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[3][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[3][i]);
    }
    
    //    4. wildtype male
    glColor3ub(150, 150, 150);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[4][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[4][i]);
    }
  
    //    5. drive male
    glColor3ub(200,   0,   0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[5][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[5][i]);
    }

    //    6. infertile male/resistant (vs other in HOMER_HS)
    glColor3ub(138,  43, 226);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[6][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[6][i]);
    }
    //    7. other male (vs rescue lost in HOMER_HS)
    glColor3ub(  0, 100,   0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[7][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[7][i]);
    }
    
    
    glEnd();
    // Outline
    glColor3ub(255, 255, 255);
    glBegin(GL_LINE_LOOP);
    glVertex2i(0,0);
    glVertex2i(Width*Scale,0);
    glVertex2i(Width*Scale,150);
    glVertex2i(0,150);
    glEnd();
    // Labels for the summary graph
    i = -40;
    /*
    glColor3ub(200, 200, 200); glRasterPos2f(i+=40, 165);  printGLf(glw->font,"f(w)");
    glColor3ub(255,   0,   0); glRasterPos2f(i+=45, 165);  printGLf(glw->font,"f(dr)");
    glColor3ub(186,  85, 211); glRasterPos2f(i+=50, 165);  printGLf(glw->font,"f(r1)");
    glColor3ub(  0, 150,   0); glRasterPos2f(i+=55, 165);  printGLf(glw->font,"f(r2)");
    glColor3ub(150, 150, 150); glRasterPos2f(i+=62, 165);  printGLf(glw->font,"m(w)");
    glColor3ub(200,   0,   0); glRasterPos2f(i+=63, 165);  printGLf(glw->font,"m(dr)");
    glColor3ub(138,  43, 226); glRasterPos2f(i+=64, 165);  printGLf(glw->font,"m(r1)");
    glColor3ub(  0, 100,   0); glRasterPos2f(i+=65, 165);  printGLf(glw->font,"m(r2)");
    */
    glColor3ub(200, 200, 200); glRasterPos2f(i+=40, 165);  printGLf(glw->font,"f(w)");     //0
    glColor3ub(255,   0,   0); glRasterPos2f(i+=45, 165);  printGLf(glw->font,"f(dr)-f");  //1
    glColor3ub(186,  85, 211); glRasterPos2f(i+=50, 165);  printGLf(glw->font,"f(r1)");    //2
    glColor3ub(  0, 150,   0); glRasterPos2f(i+=55, 165);  printGLf(glw->font,"f(r2)-if"); //3
    glColor3ub(150, 150, 150); glRasterPos2f(i+=62, 165);  printGLf(glw->font,"m(w)");     //4
    glColor3ub(200,   0,   0); glRasterPos2f(i+=63, 165);  printGLf(glw->font,"m(dr)");    //5
    glColor3ub(138,  43, 226); glRasterPos2f(i+=64, 165);  printGLf(glw->font,"m(r1)");    //6
    glColor3ub(  0, 100,   0); glRasterPos2f(i+=65, 165);  printGLf(glw->font,"r2");       //7
      
  } else {
    // Initialize buckets
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<9; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    // Fill in buckets
    //printf("%llu\n",GenMax);
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));
      //    0. wildtype female
      //    1. drive female
      //    2. semi fertile female
      //    3. drive male
      //    4. semi fertile carrier male
      //    5. wildtype male
      
      //printf("0: %llu\n",  Gen[i].geno_ze);
      //printf("1: %llu\n",  Gen[i].geno_on);
      //printf("2: %llu\n",  Gen[i].geno_tw);
      //printf("3: %llu\n",  Gen[i].geno_th);
      //printf("4: %llu\n",  Gen[i].geno_fo);
      //printf("5: %llu\n",  Gen[i].geno_fi);
      //printf("*N: %llu\n", Gen[i].Popsize);
      buckets[0][j] += ((double)(Gen[i].geno_ze))/((double)(Gen[i].Popsize));
      buckets[1][j] += ((double)(Gen[i].geno_on))/((double)(Gen[i].Popsize));
      buckets[2][j] += ((double)(Gen[i].geno_tw))/((double)(Gen[i].Popsize));
      buckets[3][j] += ((double)(Gen[i].geno_th))/((double)(Gen[i].Popsize));
      buckets[4][j] += ((double)(Gen[i].geno_fo))/((double)(Gen[i].Popsize));
      buckets[5][j] += ((double)(Gen[i].geno_fi))/((double)(Gen[i].Popsize));
      buckets[6][j] += ((double)(Gen[i].geno_si))/((double)(Gen[i].Popsize));
      buckets[7][j] += ((double)(Gen[i].geno_se))/((double)(Gen[i].Popsize));
      // Just count; used for scaling
      buckets[8][j]++;
    }
    // Scale buckets
    for(i=0; i<9; i++)
      for(j=0; j<nb; j++) {
        buckets[i][j] = 150-((buckets[i][j]/buckets[8][j])*150);
      }

    // Jump tp drawing section to cause the graph to actually be drawn
    goto draw;
  }
}

void DrawGlobalPop(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0; static double *buckets[2];

  if ( nGen <= 1 )
    return;

  // Translate
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
    // Draw graph here from buckets
    glBegin(GL_LINES);
    // fitness
    glColor3ub(150, 150,  150);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
    }
    
    glEnd();
    // Outline
    glColor3ub(255, 255, 255);
    glBegin(GL_LINE_LOOP);
    glVertex2i(0,0);
    glVertex2i(Width*Scale,0);
    glVertex2i(Width*Scale,150);
    glVertex2i(0,150);
    glEnd();
    // Labels for the summary graph
    i = -40;
    glColor3ub(150,   0, 150);  glRasterPos2f(i+=40, 165);

  } else {
    // Initialize buckets
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<2; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    // Fill in buckets
    //printf("%llu\n",GenMax);
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));
      buckets[0][j] += ((double)(Gen[i].Popsize))/((double)(GlobalMax));
      //buckets[0][j] += ((double)(Gen[i].Popsize))/((double)(GenMax));
      // Just count; used for scaling
      buckets[1][j]++;
    }
    // Scale buckets
    for(i=0; i<2; i++)
      for(j=0; j<nb; j++) {
        buckets[i][j] = 150-((buckets[i][j]/buckets[1][j])*150);
      }

    // Jump tp drawing section to cause the graph to actually be drawn
    goto draw;
  }
}


void DrawRangeMean(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0;  static double *buckets[17];

  if ( nGen <= 1 )
    return;

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
      /* Draw graph here from buckets */
      glBegin(GL_LINES);
  glColor3ub(60, 60, 60);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
  }
  glColor3ub(255, 0, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
  }
  glColor3ub(255, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
  }
  glColor3ub(0, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[3][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[3][i]);
  }
  glColor3ub(100, 0, 100);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[4][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[4][i]);
  }
  glColor3ub(0, 100, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[5][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[5][i]);
  }
  glColor3ub(75, 75, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[6][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[6][i]);
  }
  glColor3ub(0, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[7][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[7][i]);
  }
  glColor3ub(255, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[8][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[8][i]);
  }
  glColor3ub(0, 180, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[9][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[9][i]);
  }
  glColor3ub(0, 255, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[10][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[10][i]);
  }
  glColor3ub(0, 200, 55);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[11][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[11][i]);
  }
  glColor3ub(100, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[12][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[12][i]);
  }
  glColor3ub(0, 255, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[13][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[13][i]);
  }
  glColor3ub(50, 130, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[14][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[14][i]);
  }
  glColor3ub(255,255,255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[15][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[15][i]);
  }

  glEnd();
    // Outline
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(0,0);
  glVertex2i(Width*Scale,0);
  glVertex2i(Width*Scale,150);
  glVertex2i(0,150);
  glEnd();
    // Labels for the summary graph
  i = -40;

  } else {
    // Initialize buckets
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<17; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    // Fill in buckets
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));
      int l;
      for(l=0; l<((int)(pow((double)(2),(double)(nTraits)))); l++){
        // Current Fitness; red
        buckets[l][j] += (((double)(Gen[i].range.meanr[l]))/(((double)(Width))*((double)(Height))));
        //buckets[l][j] += ((double)(Gen[i].range.meanr[l]));
        //printf("\n----%lf",((double)(Gen[i].range.meanr[l]))/(((double)(Width))*((double)(Height))));
      }
      // Just count; used for scaling
      buckets[16][j]++;
    }
    // Scale buckets
    for(i=0; i<17; i++)
      for(j=0; j<nb; j++)
        buckets[i][j] = 150-((buckets[i][j]/buckets[16][j])*150);

    // Jump t0 drawing section to cause the graph to actually be drawn
    goto draw;
  }
}

void DrawRangeMax(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0;  static double *buckets[17];

  if ( nGen <= 1 )
    return;

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
      /* Draw graph here from buckets */
      glBegin(GL_LINES);
  glColor3ub(60, 60, 60);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
  }
  glColor3ub(255, 0, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
  }
  glColor3ub(255, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
  }
  glColor3ub(0, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[3][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[3][i]);
  }
  glColor3ub(100, 0, 100);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[4][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[4][i]);
  }
  glColor3ub(0, 100, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[5][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[5][i]);
  }
  glColor3ub(75, 75, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[6][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[6][i]);
  }
  glColor3ub(0, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[7][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[7][i]);
  }
  glColor3ub(255, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[8][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[8][i]);
  }
  glColor3ub(0, 180, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[9][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[9][i]);
  }
  glColor3ub(0, 255, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[10][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[10][i]);
  }
  glColor3ub(0, 200, 55);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[11][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[11][i]);
  }
  glColor3ub(100, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[12][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[12][i]);
  }
  glColor3ub(0, 255, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[13][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[13][i]);
  }
  glColor3ub(50, 130, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[14][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[14][i]);
  }
  glColor3ub(255,255,255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[15][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[15][i]);
  }

  glEnd();
    // Outline
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(0,0);
  glVertex2i(Width*Scale,0);
  glVertex2i(Width*Scale,150);
  glVertex2i(0,150);
  glEnd();
    // Labels for the summary graph
  i = -40;
  } else {
    /* Initialize buckets */
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<17; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    /* Fill in buckets */
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));

      int l;
      for(l=0; l<((int)(pow((double)(2),(double)(nTraits)))); l++){
        // Current Fitness; red
        buckets[l][j] += (((double)(Gen[i].range.mrange[l]))/(((double)(Width))*((double)(Height))));
        //buckets[l][j] += ((double)(Gen[i].range.meanr[l]));
        //printf("\n----%lf",((double)(Gen[i].range.meanr[l]))/(((double)(Width))*((double)(Height))));
      }
      buckets[16][j]++;
    }
    /* Scale buckets */
    for(i=0; i<17; i++)
      for(j=0; j<nb; j++)
        buckets[i][j] = 150-((buckets[i][j]/buckets[16][j])*150);

    /* Jump t0 drawing section to cause the graph to actually be drawn */
    goto draw;
  }
}

void DrawRangeMaxSp(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0;  static double *buckets[5];

  if ( nGen <= 1 )
    return;

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
  /* Draw graph here from buckets */
  glBegin(GL_LINES);

  glColor3ub(255, 0, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
  }
  glColor3ub(255, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
  }
  glColor3ub(100, 0, 100);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
  }
  glColor3ub(255, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[3][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[3][i]);
  }

  glEnd();
    // Outline
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(0,0);
  glVertex2i(Width*Scale,0);
  glVertex2i(Width*Scale,150);
  glVertex2i(0,150);
  glEnd();
    // Labels for the summary graph
  i = -40;

  } else {
    /* Initialize buckets */
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<5; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    /* Fill in buckets */
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));
      buckets[0][j] += ((double)(Gen[i].range.mrange[1]))/(((double)(Width))*((double)(Height)));
      buckets[1][j] += ((double)(Gen[i].range.mrange[2]))/(((double)(Width))*((double)(Height)));
      buckets[2][j] += ((double)(Gen[i].range.mrange[4]))/(((double)(Width))*((double)(Height)));
      buckets[3][j] += ((double)(Gen[i].range.mrange[8]))/(((double)(Width))*((double)(Height)));

      // Just count; used for scaling
      buckets[4][j]++;
    }
    /* Scale buckets */
    for(i=0; i<5; i++)
      for(j=0; j<nb; j++)
        buckets[i][j] = 150-((buckets[i][j]/buckets[4][j])*150);

    /* Jump t0 drawing section to cause the graph to actually be drawn */
    goto draw;
  }
}

void DrawRangeMaxGn(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0;  static double *buckets[12];

  if ( nGen <= 1 )
    return;

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
      /* Draw graph here from buckets */
      glBegin(GL_LINES);

  glColor3ub(0, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
  }
  glColor3ub(0, 100, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
  }
  glColor3ub(75, 75, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
  }
  glColor3ub(0, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[3][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[3][i]);
  }
  glColor3ub(0, 180, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[4][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[4][i]);
  }
  glColor3ub(0, 255, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[5][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[5][i]);
  }
  glColor3ub(0, 200, 55);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[6][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[6][i]);
  }
  glColor3ub(100, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[7][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[7][i]);
  }
  glColor3ub(0, 255, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[8][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[8][i]);
  }
  glColor3ub(50, 130,  0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[9][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[9][i]);
  }
  glColor3ub(255,255,255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[10][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[10][i]);
  }

  glEnd();
    // Outline
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(0,0);
  glVertex2i(Width*Scale,0);
  glVertex2i(Width*Scale,150);
  glVertex2i(0,150);
  glEnd();
    // Labels for the summary graph
  i = -40;


  } else {
    /* Initialize buckets */
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<12; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    /* Fill in buckets */
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));
      buckets[0][j] += ((double)(Gen[i].range.mrange[3]))/(((double)(Width))*((double)(Height)));
      buckets[1][j] += ((double)(Gen[i].range.mrange[5]))/(((double)(Width))*((double)(Height)));
      buckets[2][j] += ((double)(Gen[i].range.mrange[6]))/(((double)(Width))*((double)(Height)));
      buckets[3][j] += ((double)(Gen[i].range.mrange[7]))/(((double)(Width))*((double)(Height)));
      buckets[4][j] += ((double)(Gen[i].range.mrange[9]))/(((double)(Width))*((double)(Height)));
      buckets[5][j] += ((double)(Gen[i].range.mrange[10]))/(((double)(Width))*((double)(Height)));
      buckets[6][j] += ((double)(Gen[i].range.mrange[11]))/(((double)(Width))*((double)(Height)));
      buckets[7][j] += ((double)(Gen[i].range.mrange[12]))/(((double)(Width))*((double)(Height)));
      buckets[8][j] += ((double)(Gen[i].range.mrange[13]))/(((double)(Width))*((double)(Height)));
      buckets[9][j] += ((double)(Gen[i].range.mrange[14]))/(((double)(Width))*((double)(Height)));
      buckets[10][j] += ((double)(Gen[i].range.mrange[15]))/(((double)(Width))*((double)(Height)));
      
      // Just count; used for scaling
      buckets[11][j]++;
    }
    /* Scale buckets */
    for(i=0; i<12; i++)
      for(j=0; j<nb; j++)
        buckets[i][j] = 150-((buckets[i][j]/buckets[11][j])*150);

    /* Jump t0 drawing section to cause the graph to actually be drawn */
    goto draw;
  }
}
void DrawRangeNum(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0;  static double *buckets[17];

  if ( nGen <= 1 )
    return;

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
  draw:
      /* Draw graph here from buckets */
      glBegin(GL_LINES);
  glColor3ub(60, 60, 60);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
  }
  glColor3ub(255, 0, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
  }
  glColor3ub(255, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
  }
  glColor3ub(0, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[3][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[3][i]);
  }
  glColor3ub(100, 0, 100);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[4][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[4][i]);
  }
  glColor3ub(0, 100, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[5][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[5][i]);
  }
  glColor3ub(75, 75, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[6][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[6][i]);
  }
  glColor3ub(0, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[7][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[7][i]);
  }
  glColor3ub(255, 100, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[8][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[8][i]);
  }
  glColor3ub(0, 180, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[9][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[9][i]);
  }
  glColor3ub(0, 255, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[10][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[10][i]);
  }
  glColor3ub(0, 200, 55);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[11][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[11][i]);
  }
  glColor3ub(100, 0, 255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[12][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[12][i]);
  }
  glColor3ub(0, 255, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[13][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[13][i]);
  }
  glColor3ub(50, 130, 0);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[14][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[14][i]);
  }
  glColor3ub(255,255,255);
  for(i=1; i<nb; i++) {
    glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[15][i-1]);
    glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[15][i]);
  }

  glEnd();
    // Outline
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(0,0);
  glVertex2i(Width*Scale,0);
  glVertex2i(Width*Scale,150);
  glVertex2i(0,150);
  glEnd();
    // Labels for the summary graph
  i = -40;

  } else {
    /* Initialize buckets */
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<17; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    /* Fill in buckets */
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));

      int l;
      for(l=0; l<((int)(pow((double)(2),(double)(nTraits)))); l++){
        // Current Fitness; red
        buckets[l][j] += 10*(((double)(Gen[i].range.noisolate[l]))/(((double)(Width))*((double)(Height))));
        //buckets[l][j] += ((double)(Gen[i].range.meanr[l]));
        //printf("\n----%lf",((double)(Gen[i].range.meanr[l]))/(((double)(Width))*((double)(Height))));
      }

      // Just count; used for scaling
      buckets[16][j]++;
    }
    /* Scale buckets */
    for(i=0; i<17; i++)
      for(j=0; j<nb; j++)
        buckets[i][j] = 150-((buckets[i][j]/buckets[16][j])*150);

    /* Jump t0 drawing section to cause the graph to actually be drawn */
    goto draw;
  }
}

void DrawPopSize(GLWindow *glw, float tx, float ty)
{
  int i,j;  static int nb=0;  static double *buckets[17];

  if ( nGen <= 1 )
    return;

  /* Translate */
  glLoadIdentity();
  glTranslatef(tx, ty, 1.0f);

  if(nb) {
    draw:
    /* Draw graph here from buckets */
    glBegin(GL_LINES);
    glColor3ub(60, 60, 60);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[0][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[0][i]);
    }
    glColor3ub(255, 0, 0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[1][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[1][i]);
    }
    glColor3ub(255, 0, 255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[2][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[2][i]);
    }
    glColor3ub(0, 0, 255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[3][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[3][i]);
    }
    glColor3ub(100, 0, 100);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[4][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[4][i]);
    }
    glColor3ub(0, 100, 255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[5][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[5][i]);
    }
    glColor3ub(75, 75, 255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[6][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[6][i]);
    }
    glColor3ub(0, 100, 0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[7][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[7][i]);
    }
    glColor3ub(255, 100, 0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[8][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[8][i]);
    }
    glColor3ub(0, 180, 255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[9][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[9][i]);
    }
    glColor3ub(0, 255, 255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[10][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[10][i]);
    }
    glColor3ub(0, 200, 55);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[11][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[11][i]);
    }
    glColor3ub(100, 0, 255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[12][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[12][i]);
    }
    glColor3ub(0, 255, 0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[13][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[13][i]);
    }
    glColor3ub(50, 130, 0);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[14][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[14][i]);
    }
    glColor3ub(255,255,255);
    for(i=1; i<nb; i++) {
      glVertex2d((((double)(i-1))/(nb-1))*Width*Scale, buckets[15][i-1]);
      glVertex2d((((double)(i))/(nb-1))*Width*Scale,   buckets[15][i]);
    }

    glEnd();
    // Outline
    glColor3ub(255, 255, 255);
    glBegin(GL_LINE_LOOP);
    glVertex2i(0,0);
    glVertex2i(Width*Scale,0);
    glVertex2i(Width*Scale,150);
    glVertex2i(0,150);
    glEnd();
    // Labels for the summary graph
    i = -40;

  } else {
    /* Initialize buckets */
    nb = ((200<nGen)?(200):(nGen));
    for(i=0; i<17; i++) {
      buckets[i] = malloc(nb*sizeof(double));
      memset(buckets[i], 0, nb*sizeof(double));
    }
    /* Fill in buckets */
    for(i=0; i<nGen; i++) {
      j = (int)((float)(((((double)(i))/(nGen-1))*(nb-1))));
      int l;
      for(l=0; l<((int)(pow((double)(2),(double)(nTraits)))); l++){
        // Current Fitness; red
        buckets[l][j] += (((double)(Gen[i].range.sppop[l]))/(Max*((double)(Width))*((double)(Height))));
        //buckets[l][j] += ((double)(Gen[i].range.meanr[l]));
        //printf("\n----%lu",((double)(Gen[i].range.sppop[l])));
      }
      buckets[16][j]++;
    }
    /* Scale buckets */
    for(i=0; i<17; i++)
      for(j=0; j<nb; j++)
        buckets[i][j] = 150-((buckets[i][j]/buckets[16][j])*150);

    /* Jump t0 drawing section to cause the graph to actually be drawn */
    goto draw;
  }
}

void DrawOutlineTranslucentBox(int x1, int y1, int x2, int y2)
{
  // Translucent background
  glColor4ub(0, 0, 0, 128);
  glBegin(GL_QUADS);
  glVertex2i(x1,y1);
  glVertex2i(x2,y1);
  glVertex2i(x2,y2);
  glVertex2i(x1,y2);
  glEnd();
  // Outline
  glColor3ub(255, 255, 255);
  glBegin(GL_LINE_LOOP);
  glVertex2i(x1,y1);
  glVertex2i(x2,y1);
  glVertex2i(x2,y2);
  glVertex2i(x1,y2);
  glEnd();
  glEndList();
}

void DrawHoverBox(GLWindow *glw)
{
  int i,j,k,x,y;  char b[1024],buf[4096];

  /* Just return if out of bounds */
  if( ((i=(x=Mx)/((int)Scale)) >= Width) || ((j=(y=My)/((int)Scale)) >= Height) )
    return;

  // Draw the hover box background
  if( Gen[Frame-1].p[i][j].ni ) {
    // Bound and translate
    if( (y+150) > (Height*Scale) ) y = (Height*Scale)-150;
    if( (x+200) > (Width*Scale)  ) x = (Width*Scale) -200;
    glLoadIdentity();
    glTranslatef(x, y, 1.0f);
#if PRESENTATION_MODE
    DrawOutlineTranslucentBox(0,0,250,150);
#else
    DrawOutlineTranslucentBox(0,0,150,60);
#endif
  } else {
  // Bound and translate
    if( (y+50)  > (Height*Scale) ) y = (Height*Scale)-50;
    if( (x+150) > (Width*Scale)  ) x = (Width*Scale) -150;
    glLoadIdentity();
    glTranslatef(x, y, 1.0f);
#if PRESENTATION_MODE
    DrawOutlineTranslucentBox(0,0,200,100);
#else
    DrawOutlineTranslucentBox(0,0,150,50);
#endif
  }

  /* Less info for empty patches */
  if( !Gen[Frame-1].p[i][j].ni ) {
    glColor3ub(255,255,255);
    glRasterPos2f(15.0f,15.0f);
    printGLf(glw->font,"Patch popsize: 0");
    //glRasterPos2f(15.0f,35.0f);
    //strcpy(buf, "env: ");
    //for(k=0; k<nTraits; k++) {
      //sprintf(b, "%d ", Gen[Frame-1].p[i][j].l[k]);
      //strcat(buf,b);
    //}
    //printGLf(glw->font,"%s",buf);
  } else {
    /* Normal Stats */
    glColor3ub(255,255,255);
    glRasterPos2f(15.0f,15.0f);
    printGLf(glw->font,"Patch popsize: %llu", Gen[Frame-1].p[i][j].ni);
#if HOMER || HOMER_HS || HOMER_KO
    glRasterPos2f(15.0f,35.0f);
    printGLf(glw->font,"Gene-drive: %.3lf", Gen[Frame-1].p[i][j].f);
    //glRasterPos2f(15.0f,55.0f);
    //printGLf(glw->font,"Prolactin: %.3lf", Gen[Frame-1].p[i][j].f);
#endif
#if HOMINGDRIVE
    glRasterPos2f(15.0f,35.0f);
    //printGLf(glw->font,"Prolactin: %.3lf", Gen[Frame-1].p[i][j].f);
#endif

    
  }

}

/***************************************************************************************
 * Scene drawing code.  Component layout mostly
 ***************************************************************************************/

static void DrawScene(GLWindow *glw)
{
  Plot.x = 1;
  Plot.y = 480;
  Plot.w = 500;
  Plot.h = 200;

  /* Ugly hackish bugfix */
  if(!Frame) Frame = 1;
  if(Frame > nGen) Frame = nGen;

  /* Clear the old scene */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Attache 3d plot draw function here
  //Draw3DDeme                    (glw);

  /* Draw the deme graph */
  DrawDemeGraph                 (glw, 0.5f, 1.0f);


  DrawSummaryGraph              (glw, Width*Scale+20.0f, 10.0f);
  DrawGenotypePop               (glw, Width*Scale+20.0f, 180.0f);
  DrawGlobalPop                 (glw, Width*Scale+20.0f, 350.0f);
  DrawGenerationSlider          (glw, Width*Scale+20.0f, 510.0f);
  DrawGenotypeHistogram         (glw, Width*Scale+20.0f, 560.0f);
    
#if READ_DISTANCE
  DrawDVD                       (glw, 210.0f+NICHES*(20+HS), 420.0f+NICHES*(20+HS));
#endif

  // Draw the text-based statistics
  DrawTextStats                 (glw, 0.5f, Height*Scale+40.0f);
  //DrawTextStats               (glw, Width*Scale+20.0f, 340.0f);
  //DrawSpeciesLegend           (glw, 0.0f, Height*Scale+220.0f);
  //DrawSpeciesLegend           (glw, 0.0f, Height*Scale+50.0f);

  // Draw gneration slider
  DrawGenerationSlider          (glw, 0.5f, Height*Scale+11.0f);
  //DrawGenerationSlider        (glw, Width*Scale+20.0f, Height*Scale+340.0f);
    

#if !SIMPLE
  /* Draw hover-box */
  DrawHoverBox(glw);

  /* Draw summary graph */
  //DrawSummaryGraph            (glw, 0.0f, Height*Scale+45.0f);
  //DrawRangeMean               (glw, Width*Scale+20.0f, 1.0f);
  //DrawMatch                   (glw, Width*Scale+20.0f, 1.0f);
  //DrawFitness                 (glw, Width*Scale+20.0f, 1.0f);
  //DrawRangeMax                  (glw, Width*Scale+20.0f, 1.0f);
  //DrawRangeMaxSp              (glw, Width*Scale+20.0f, 310.0f);
  //DrawRangeMaxGn              (glw, Width*Scale+20.0f, 465.0f);
  //DrawRangeNum                (glw, Width*Scale+20.0f, 620.0f);
  //DrawPopSize                 (glw, Width*Scale+20.0f, 155.0f);
  //DrawRangeNum                (glw, Width*Scale+20.0f, 310.0f);
#endif

  /* Draw the niche histogram */

    /* Draw preference graph (Cramer's graph) */
  //  DrawPreferenceGraph           (glw, Width*Scale+500.0f, 210.0f);
  //  DrawNicheHistogram            (glw, Width*Scale+780.0f, 200.0f+NICHES*(20+HS) + 20);
  
  /* Swap to display everything */
  glXSwapBuffers(glw->dpy, glw->win);
}

/***************************************************************************************
 * Window creation and event handler
 ***************************************************************************************/

void Deme_Down(const int x, const int y)
{
  // We only care about mouse down events within the widget area
    if( (x > 1) && (x < 501) && 
        (y > 480) && (y < 980) ) {
        // Mouse rotation controll
        if( !Plot.d ) {
          Plot.d  = 1;
          Plot.ox = Plot.nx = x;
          Plot.oy = Plot.ny = y;
        }
  }
}

void Deme_Up(const int x, const int y)
{
  if( Plot.d ) {
    Plot.d  = 0;
    Plot.nx = x;
    Plot.ny = y;
    Plot.rot.z += Plot.nx - Plot.ox;
    Plot.rot.x += Plot.ny - Plot.oy;
    Plot.ox = Plot.nx;
    Plot.oy = Plot.ny;
  }
}

void Deme_Move(const int x, const int y)
{
  // Mouse rotation controll
  if( Plot.d ) {
    Plot.nx = x;
    Plot.ny = y;
  }
}

void Graph()
{
  XEvent event;  GLWindow glw;  KeySym key;  
  char keys[255],d=0,b[1024],b_nGen[1024];  char *title = "niche_view";
  double st=0,ct=0;  struct timeval tv;
  int i,t,done=0,x=0,y=0,in=0,oin=0,update=1;
  FILE *gnuplot;

  /* Create window and initialize it */
  Frame  = (int)nGen;
  glw.id = 0;
  CreateGLWindow(&glw, title, X, Y);
  InitGLWindow(&glw);

  /* Process events */
  while (!done){
    if (XCheckIfEvent(glw.dpy, &event, AnyEvent, NULL)){
      switch (event.type) {
      case Expose:
	if (event.xexpose.count != 0)
	  break;
	update = 1;
	break;
      case ConfigureNotify:
	if ((event.xconfigure.width != glw.width) || (event.xconfigure.height != glw.height)) {
	  glw.width = event.xconfigure.width;
	  glw.height = event.xconfigure.height;
	  ResizeGLScene(event.xconfigure.width, event.xconfigure.height);
	}
	break;
      case ClientMessage:    
	if (*XGetAtomName(glw.dpy, event.xclient.message_type) == *"WM_PROTOCOLS")
	  done = 1;
	break;
      case KeyPress:
	if( XLookupString(&event.xkey, keys, 255, &key, 0) == 1 ) {
	  switch (*keys){
	  case '+': 
	    Speed++; 
	    if (Speed <  0) Speed =  0;
	    if (Speed > 60) Speed = 60;
	    update = 1;
	    break;
	  case '-':
	    Speed--; 
	    if (Speed <  0) Speed =  0;
	    if (Speed > 60) Speed = 60;
	    update = 1;
	    break;
	  case '<':
	    if( --Scale < 1 )
	      Scale = 1;
	    update = 1;
	    break;
	  case '>':
	    Scale++;
	    update = 1;
	    break;
	  case ',':
	    if(Frame > 1) Frame--;
	    update = 1;
	    break;
	  case '.':
	    if(Frame < nGen) Frame++;
	    update = 1;
	    break;
	  case ' ':
	    Speed = 0;
	    update = 1;
	    break;
	  case 'c':
	    ColorMode ^= 1;
	    update = 1;
	    break;
	  }
	}
	break;
      case ButtonRelease:
	// Attache 3d plot button here
	Deme_Up(event.xbutton.x,event.xbutton.y);
	d = 0;
	break;
      case MotionNotify:
	// Attache 3d plot button here
	Deme_Move(event.xbutton.x,event.xbutton.y);
	Mx = event.xbutton.x;
	My = event.xbutton.y;
	if(d) {
	  Trans[0] += event.xbutton.x-x;
	  Trans[1] += event.xbutton.y-y;
	  x = event.xbutton.x;
	  y = event.xbutton.y;
	  Frame = (((double)x)/(Width*Scale))*nGen;
	  update = 1;
	}
	if( (Mx >= 0) && (Mx < Width*Scale) && (My >=0) && (My < Height*Scale) ) {
	  in = 1;
	  update = 1;
	} else {
	  in = 0;
	}
	if( oin != in ) {
	  update = 1;
	  oin = in;
	}
	break;
      case ButtonPress:
	// Attache 3d plot button here
	Deme_Down(event.xbutton.x,event.xbutton.y);
	switch (event.xbutton.button){
	case 1:
	  x = event.xbutton.x; 
	  y = event.xbutton.y;
	  // Mouse down over the niche histogram has special meaning
	  if( (x > Width*Scale+20.0f) && (x < NICHES*(20+HS)+Width*Scale+20.0f) &&
	      (y > 220) && (y<320) ) {
	    Filter[(int)((x-(Width*Scale+20.0f))/(20+HS))] ^= 1;
	    update = 1;
	  }
	  // Slider
	  if( (x < Width*Scale) && (y > Height*Scale+10) && (y < Height*Scale+40) ) {
	    d = 1;
	    Frame = (((double)x)/(Width*Scale))*nGen;
	    update = 1;
	  }
	  break;
	case 3: 
	  Frame = 0;
	  gettimeofday(&tv,NULL);
	  st = tv.tv_sec + tv.tv_usec/1000000.;
	  break;
	case 4:
	  if( --Scale < 1 )
	    Scale = 1;
	  update = 1;
	  break;
	case 5:
	  Scale++;
	  update = 1;
	  break;
	}
      }
    } else {
      /* Dump if needed and exit */
      if(dump) {
	if(dump < 0) dump = nGen-1;
	t = Frame;
	for(i=0; i<nGen; i++) {
	  if( !dump || ( !(i%dump) && (i || dump!=nGen-1) ) ) {
	    Frame = i+1;
	    if( !(PSFile0 = PS_Start("spc.eps",(nTraits)/2.+.25,(nTraits)/2.+.25)) ) {
	      fprintf(stderr,"Could not open output post-script file! (tmp1.eps)\n");
	      exit(1);
	    }
	    if( !(PSFile1 = PS_Start("deme.eps",Width/2.,Height/2.)) ) {
	      fprintf(stderr,"Could not open output post-script file! (tmp0.eps)\n");
	      exit(1);
	    }
	    DrawScene(&glw);
	    sprintf(b,"%llu",Gen[Frame-1].g);
	    sprintf(b_nGen,"%llu",nGen);
	    t = strlen(b);
	    b[0] = 0;
	    while(t++<strlen(b_nGen))
	      strcat(b,"0");
	    sprintf(b_nGen,"%s",b);
	    PS_End(PSFile0);
	    PS_End(PSFile1);
	    PSFile0 = NULL;
	    PSFile1 = NULL;
	  }
	}
	Frame = t;
	exit(0);
      }
      if ( (Frame < nGen) && (Speed) ){
	gettimeofday(&tv,NULL);
	ct = tv.tv_sec + tv.tv_usec/1000000.;
	if( (!d) && (ct-st > 1./Speed) ) {
	  Frame++;
	  DrawScene(&glw);
	  update = 0;
	  st = ct;
	} else {
	  /* Not ready to draw next frame yet / dragging */
	  // Updated if needed, otherwise yeild.
	  if( update ) {
	    DrawScene(&glw);
	    update = 0;
	  } else {
	    DrawScene(&glw);
	    yeild();
	  }
	}
      } else {
	/* Not ready to draw next frame yet / stopped */
	// Updated if needed, otherwise yeild.
	if( update ) {
	  DrawScene(&glw);
	  update = 0;
	} else {
	  DrawScene(&glw);
	  yeild();
	}
      }
    }

  }

  // We are done, kill the window
  KillGLWindow(&glw);
}

/***************************************************************************************
 * Application entry point
 ***************************************************************************************/

void Dumpall()
{
  if( nGen ) {
    Frame = nGen;
  } else {
    printf("error: frame < 0.\n");  
    return;
  }

  printf("Cycle:%llui\n",           Gen[Frame-1].g);        // "Mating cycle" (old version: Generation)
  printf("pop:%llui\n",             Gen[Frame-1].Popsize);  // PopSize
  printf("Tcolonize:%di\n",         Timetocolonize);        // colonization time (in cycles)
  printf("Nmax:%di\n",              GlobalMax);             // MaxPopSize
  printf("Nmin:%di\n",              GlobalMin);             // MinPopSize
  printf("Tmin:%di\n",              TMinPop);               // T for MinPopSize
  printf("Treduction:%di\n",        TPopSizeReduction);     // time when popsize is <0.01(GlobalMax) (in cycles)
  WriteDemes();
}

void UseageError()
{
  fprintf(stderr,"Useage\n\tview <dataset> [-x <width>] [-y <height>] [-dump [interval]]\n");
  exit(1);
}

int main(int argc, char **argv)
{
  int i;

  /* Check command line args and read data file */
  if ( (argc < 2) || (argc > 7))
    UseageError();
  FNPrefix = strdup(argv[1]);
  for(i=2; i<argc; i++) {
    if(!strcmp(argv[i],"-x")) {
      /* width of window */
      if(i+1 >= argc)
        UseageError();
      if(sscanf(argv[i+1],"%d",&X) != 1)
        UseageError();
      i++;
    }
    if(!strcmp(argv[i],"-y")) {
      /* height of window */
      if(i+1 >= argc)
        UseageError();
      if(sscanf(argv[i+1],"%d",&Y) != 1)
        UseageError();
      i++;
    }
    if(!strcmp(argv[i],"-dump")) {
      dump = -1;
      /* silent/dump mode */
      if(i+1 < argc)
        if(sscanf(argv[i+1],"%d",&dump) == 1)
          i++;
    }
    if(!strcmp(argv[i],"-amc")) {
      // Dump the number of generations until amc >= .5
      if( (i=readDataFile(argv[1])) ) {
        printf("-1\n");
        exit(i);
      }
      printf("%d\n",Amcg);
      return 0;
    }
    if(!strcmp(argv[i],"-gens")) {
      // Dump the number of generations in the data file
      if( (i=readDataFile(argv[1])) ) {
        printf("-1\n");
        exit(i);
      }
      printf("%llu\n",Gen[nGen-1].g);
      return 0;
    }
    if(!strcmp(argv[i],"-clean")) {
      // Dump a flag for whether or not the data files represent a run that exited cleanly.
      // Really, this dumps 1 if the trailer is present or 0 otherwise
      if( (i=readDataFile(argv[1])) ) {
        printf("-1\n");
        exit(i);
      }
      printf("%d\n",Clean);
      return 0;
    }
    if(!strcmp(argv[i],"-spcs")) {
      // Dumps the number of species in the last gen
      if( (i=readDataFile(argv[1])) ) {
        printf("-1\n");
        exit(i);
      }
      printf("%d\n",Gen[nGen-1].tns);
      return 0;
    }
    if(!strcmp(argv[i],"-fn")) {
      // Dumps the number of "filled" niches in the last gen
      if( (i=readDataFile(argv[1])) ) {
        printf("-1\n");
        exit(i);
      }
      printf("%d\n",Gen[nGen-1].fn);
      return 0;
    }
    if(!strcmp(argv[i],"-kpnvar")) {
      // Dumps the average variance of the even ... read the changelog/code ...
      if( (i=readDataFile(argv[1])) ) {
        printf("-1\n");
        exit(i);
      }
      printf("%lf\n",Gen[nGen-1].vark);
      return 0;
    }
    if(!strcmp(argv[i],"-spcspfn")) {
      // Dumps the number of species per filled niche in the last gen
      if( (i=readDataFile(argv[1])) ) {
        printf("-1\n");
        exit(i);
      }
      printf("%lf\n",Gen[nGen-1].spfn);
      return 0;
    }
    if(!strcmp(argv[i],"-dumpall")) {
      // Dumps a large number of statistics in html format
      if( (i=readDataFile(argv[1])) ) {
        printf("Read of data file failed.\n");
        exit(i);
      }

      RecordInvasion();
      RecordReduction();
      Dumpall();

      return 0;
    }
  }

  // read data
  if( (i=readDataFile(argv[1])) )
    exit(i);

  RecordInvasion();
  RecordReduction();
  
  // Run the graph
  Graph();

  // Return success
  return 0;
}
