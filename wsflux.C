/*
  Program: wsflux.C
  Author:  D. Trinkle
  Date:    2010 August 1
  Purpose: Reads in a CHGCAR (header) file, and determines the set of
           neighboring directions to be considered, along with the 
	   Wigner-Seitz facet area divided by vector length.  These prefactors
	   are then multiplied by delta rho to determine the transition matrix
	   element in the gradient flow algorithm.

  Param.:  chgcar

  Flags:   MEMORY:  not used
           VERBOSE: not used
           TESTING: usual screen diahrea

  Algo.:   After reading in the chgcar file to get the grid basis vectors,
           we construct the Wigner-Seitz cell.

  Output:  Whatever I need to output.

*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include "matrix.H"
#include "chgcar.H"

//****************************** STRUCTURES ****************************

//***************************** SUBROUTINES ****************************
inline FILE* myopenr(char *filename) 
{ if (filename[0] == '-') return stdin;
  else                    return fopen(filename, "r"); }

inline FILE* myopenw(char *filename) 
{ if (filename[0] == '-') return stdout;
  else                    return fopen(filename, "w"); }
  
inline void myclose(FILE *f) 
{ if ((f != stdin) && (f != stdout) && (f != NULL)) fclose(f); }

inline int is_blank (char c) {return (c == ' ') || (c == '\t') || (c == '\n');}

inline void nextnonblank (char* &p) 
{ for ( ; ((*p)!='\0') && is_blank(*p); ++p) ;}

inline void nextblank (char* &p) 
{ for ( ; ((*p)!='\0') && (!is_blank(*p)); ++p) ;}

const char COMMENT_CHAR = '#';
const char EOF_CHAR = '&';

inline void nextnoncomment (FILE* infile, char* dump, const int &size) 
{ do {fgets(dump, size, infile);} 
  while ((!feof(infile)) && (dump[0] == COMMENT_CHAR)); }


const double TOLER=1e-8;

// is a vector r inside of the WS cell defined by the set of vectors R?
// NOTE: R[n][0..2] = cartesian coord, R[n][3] = (R.R)/2
inline int incell (double r[3], int Nneigh, double** R) 
{
  for (int n=0; n<Nneigh; ++n) 
    if (dot(r, R[n]) > (R[n][3]+TOLER)) return 0;
  return 1;
}

// want this list in ascending order 
int vert_comp(const void* a, const void* b) 
{
  double *av=*((double**)a);
  double *bv=*((double**)b);
  if ( av[3] < bv[3] ) return -1;
  if ( av[3] > bv[3] ) return  1;
  return 0;
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "[-m MEMORY] [-hvt] chgcar";

const char* ARGEXPL = 
"  chgcar:  CHGCAR file to be used for Bader partitioning\n\
\n\
  -m [MEMORY]  suggest necessary memory for input\n\
  -h  HELP\n\
  -v  VERBOSE\n\
  -t  TESTING";

const char* FILEEXPL =
"Expects VASP CHGCAR format; doesn't read past the grid density";

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 128;  // Our default storage (not used)
  
  char ch;
  while ((ch = getopt(argc, argv, "hvtm:")) != -1) {
    switch (ch) {
    case 'm':
      MEMORY = (int)strtol(optarg, (char**)NULL, 10);
      break;
    case 'v':
      VERBOSE = 1;
      break;
    case 't':
      TESTING = 1;
      VERBOSE = 1;
      break;
    case 'h':
    case '?':
    default:
      ERROR = 1;
    }
  }
  argc -= optind;
  argv += optind;
  if (MEMORY < 1) MEMORY = 1;

  if (argc < NUMARGS) ERROR = 2;
  // All hell broken loose yet?
  if (ERROR != 0) {
    fprintf(stderr, "%s %s\n", progname, ARGLIST);
    fprintf(stderr, "%s\n", ARGEXPL);
    fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
    exit(ERROR);
  }

  // ****************************** INPUT ****************************
  // Command line parameters:
  char* infile_name = argv[0];

  // check input:
  FILE* infile = myopenr(infile_name);
  if (infile == NULL) {
    fprintf(stderr, "Could not open %s\n", infile_name);
    ERROR = 2;
  }

  // If we've thrown an error already, then get out.
  if (ERROR) exit(ERROR);

  // parse that file...
  char dump[512];
  double a0, alatt[9];
  int Ngrid[3];

  read_CHGCAR_header(infile, dump, a0, alatt);
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%d %d %d", Ngrid, Ngrid+1, Ngrid+2);

  if (TESTING) {
    printf("## a0= %.12g\n", a0);
    for (int i=0; i<3; ++i) {
      printf("## a%d=", i+1);
      for (int d=0; d<3; ++d)
	printf(" %.12g", alatt[3*d+i]);
      printf("\n");
    }
    printf("## Ngrid= %d %d %d\n", Ngrid[0], Ngrid[1], Ngrid[2]);
  }

  for (int i=0; i<3; ++i)
    for (int d=0; d<3; ++d)
      alatt[3*d+i] *= a0/(double)Ngrid[i];

  if (TESTING) {
    for (int i=0; i<3; ++i) {
      printf("## A%d=", i+1);
      for (int d=0; d<3; ++d)
	printf(" %.12g", alatt[3*d+i]);
      printf("\n");
    }
  }

  myclose(infile);
  
  // ***************************** ANALYSIS **************************
  // 0. Generate (pruned) list of neighboring vectors that bound the
  // Wigner-Seitz cell.
  // Note for future: should precompute this by making a sphere of radius
  // with the largest length vector multiplied by, say, 2.
  const int Nrange=3;	// (generally) overkill; need a way to compute this
  int Nneigh = (2*Nrange+1)*(2*Nrange+1)*(2*Nrange+1)-1;
  double **R = new double*[Nneigh];
  int **nvect = new int*[Nneigh];
  Nneigh=0;
  int nv[3];
  for (nv[0]=-Nrange; nv[0]<=Nrange; ++nv[0])
    for (nv[1]=-Nrange; nv[1]<=Nrange; ++nv[1])
      for (nv[2]=-Nrange; nv[2]<=Nrange; ++nv[2])
	if (! zero_vect(nv) ) {
	  nvect[Nneigh]=new int[3];
	  R[Nneigh]=new double[4];
	  for (int d=0; d<3; ++d) nvect[Nneigh][d]=nv[d];
	  mult_vect(alatt, nv, R[Nneigh]);
	  R[Nneigh][3] = 0.5*dot(R[Nneigh], R[Nneigh]);
	  Nneigh++;
	}
  // prune that list: basically, if R/2 isn't inside the WS cell, then
  // R isn't involved in bounding the WS cell.
  for (int n=(Nneigh-1); n>=0; --n) {
    double r[3] = {0.5*R[n][0], 0.5*R[n][1], 0.5*R[n][2]};
    if (! incell(r, Nneigh, R) ) {
      // need to remove this entry...
      delete[] nvect[n];
      delete[] R[n];
      for (int i=n; i<(Nneigh-1); ++i) {
	nvect[i] = nvect[i+1];
	R[i] = R[i+1];
      }
      --Nneigh;
    }
  }
  
  if (TESTING) {
    printf("## Nneigh= %d\n", Nneigh);
    for (int n=0; n<Nneigh; ++n) {
      printf("## %3d%3d%3d  %.12g %.12g %.12g  %.12g\n",
	     nvect[n][0], nvect[n][1], nvect[n][2],
	     R[n][0], R[n][1], R[n][2],
	     R[n][3]);
    }
  }

  // 1. Run over each facet to find all of the vertex points
  int MAXVERT = (Nneigh-2)*(Nneigh-4);
  double** rvert = new double*[MAXVERT];
  for (int nv=0; nv<MAXVERT; ++nv) rvert[nv] = new double[4]; // [3] = phi
  double alpha[Nneigh];
  int Nnonzero=Nneigh;
  for (int n=0; n<Nneigh; ++n) {
    int nvert = 0;
    double Rdot[9], detR, Rinv[9], R2[3];
    for (int d=0; d<3; ++d) Rdot[d] = R[n][d];
    R2[0] = R[n][3];
    for (int nA=0; nA<Nneigh; ++nA) {
      for (int d=0; d<3; ++d) Rdot[3+d] = R[nA][d];
      R2[1] = R[nA][3];
      for (int nB=(nA+1); nB<Nneigh; ++nB) {
	for (int d=0; d<3; ++d) Rdot[6+d] = R[nB][d];
	R2[2] = R[nB][3];
	detR = inverse(Rdot, Rinv);
	if (fabs(detR)>TOLER) {
	  mult_vect(Rinv, R2, rvert[nvert]);
	  for (int d=0; d<3; ++d) rvert[nvert][d] *= 1./detR;
	  if (incell(rvert[nvert], Nneigh, R) ) ++nvert;
	}
      }
    }
    // check to make sure none of the vertices correspond to R/2:
    int zeroarea=0;
    for (int nv=0; nv<nvert; ++nv)
      if (fabs(dot(rvert[nv], rvert[nv])-0.5*R[n][3])<TOLER) zeroarea=1;
    if (zeroarea || (nvert==0)) {
      alpha[n]=0;
      Nnonzero--;
      continue;
    }
    // Now we have a list of all the vertices for the polygon
    // defining the facet along the direction R[n].
    // Last step is to sort the list in terms of a winding angle around
    // R[n].  To do that, we define rx and ry which are perpendicular
    // to R[n], normalized, and right-handed: ry = R x rx, so that
    // rx x ry points along R[n].
    double rx[3], ry[3];
    for (int d=0; d<3; ++d) rx[d] = rvert[0][d];
    double rdRn = dot(rx, R[n])/dot(R[n],R[n]);
    for (int d=0; d<3; ++d) rx[d] -= rdRn*R[n][d];
    rdRn = sqrt(dot(rx, rx));
    for (int d=0; d<3; ++d) rx[d] *= 1./rdRn;
    crossprod(R[n], rx, ry);
    rdRn = sqrt(dot(ry, ry));
    for (int d=0; d<3; ++d) ry[d] *= 1./rdRn;
    // now compute winding angle phi_n
    for (int nv=0; nv<nvert; ++nv)
      rvert[nv][3] = atan2(dot(rvert[nv], ry), dot(rvert[nv], rx));
    // sort that list
    qsort(rvert, nvert, sizeof(double*), vert_comp);
    alpha[n] = 0;
    for (int nv=0; nv<nvert; ++nv)
      alpha[n] += tripleprod(rvert[nv], rvert[(nv+1)%nvert], R[n]);
    alpha[n] *= 0.25/R[n][3];
    if (fabs(alpha[n])<TOLER) {
      alpha[n]=0;
      --Nnonzero;
    }
  }
  // garbage collection
  for (int nv=0; nv<MAXVERT; ++nv) delete[] rvert[nv];
  delete[] rvert;

  // ****************************** OUTPUT ***************************
  printf("%d\t#Nneigh\n", Nnonzero);
  for (int n=0; n<Nneigh; ++n)
    if (alpha[n]!=0) {
      printf("%3d%3d%3d %.15le",
	     nvect[n][0], nvect[n][1], nvect[n][2],
	     alpha[n]);
      if (TESTING)
	printf(" # %.8g %.8g %.8g\n", R[n][0], R[n][1], R[n][2]);
      else
	printf("\n");
    }
  
  // ************************* GARBAGE COLLECTION ********************
  for (int n=0; n<Nneigh; ++n) {
    delete[] nvect[n];
    delete[] R[n];
  }
  delete[] nvect;
  delete[] R;
  
  return 0;
}
