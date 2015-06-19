/*
  Program: weight_int.C
  Author:  D. Trinkle
  Date:    2010 December 28
  Purpose: Reads in a CHGCAR file, and determines the basins of attraction (either
           tied to atoms or not--set with a switch), and then integrates over
	   the basins using the weight method.  Trying to do this in the most
	   computationally efficient way possible.

	   Implementation of algorithm from:
	   “Accurate and efficient algorithm for Bader charge integration”
	   M. Yu and D. R. Trinkle. J. Chem. Phys. 134, 064111 (2011)
	   doi:10.1063/1.3553716

  Param.:  chgcar [grid1] [grid2] ...

  Flags:   MEMORY:  not used
           VERBOSE: not used
           TESTING: usual screen diahrea

  Algo.:   After reading in the chgcar file to get the grid basis vectors,
           we construct the Wigner-Seitz cell.

  Output:  Integrals over basins for each of the files; can also output 
	   weights on a similar CHGCAR-style grid.

*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include "matrix.H"
#include "chgcar.H"
#include "ws_voronoi.H"

// uncomment line below to "approximate" near grid... ONLY for testing purposes
// #define NEAR_GRID
// uncomment line below to not store the grids in sorted order; doesn't seem
// to make a big difference in runtime, though
// #define NOSORT_STORAGE

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

// convert from 3D to grid, and back
inline void tovect(const int& i, int Ngrid[3], int iv[3]) 
{
  iv[0] = i%(Ngrid[0]);
  iv[1] = (i/Ngrid[0])%(Ngrid[1]);
  iv[2] = (i/(Ngrid[0]*Ngrid[1]))%(Ngrid[2]);
}

// looks funny because we need to be able to handle PBC
inline int toindex(int Ngrid[3], int iv[3]) 
{
  return (iv[0]+Ngrid[0])%Ngrid[0] 
    + Ngrid[0]*(
		(iv[1]+Ngrid[1])%Ngrid[1]
		+ Ngrid[1]*(
			    (iv[2]+Ngrid[2])%Ngrid[2])
		);
}

// want this list in descending order--but just sort the *indices*
double* DENSITY_COMP;
// int density_comp(void* thunk, const void* a, const void* b) 
int density_comp(const void* a, const void* b) 
{
  int av=*((int*)a);
  int bv=*((int*)b);
  //  double* rho=(double*)thunk;
  double* rho=DENSITY_COMP;
  if ( rho[av] > rho[bv] ) return -1;
  if ( rho[av] < rho[bv] ) return  1;
  return 0;
}

int match_atom(const int& Natoms, double** uatom, 
	       double metric[9], double uvect[3],
	       double& r2);

inline int match_atom(const int& Natoms, double** uatom,
		      double metric[9], double uvect[3]) 
{
  double r2;
  return match_atom(Natoms, uatom, metric, uvect, r2);
}

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "[-m MEMORY] [-hvt] chgcar [grid1] [grid2] ...";

const char* ARGEXPL = 
"  chgcar:  CHGCAR file to be used for Bader partitioning\n\
  gridN:   additional grid files to integrate\n\
  -V       find Voronoi volumes instead (ignores the grid data in chgcar)\n\
  -a       do not assign basins to atoms, but give unique tags\n\
  -s       scale by total volume of cell (needed for Bader charges)\n\
  -o base  output weights to base0001\n\
  -n N     just output weight from basin N\n\
\n\
  -m [MEMORY]  suggest necessary memory for input\n\
  -h  HELP\n\
  -v  VERBOSE\n\
  -t  TESTING";

const char* FILEEXPL =
"Expects VASP CHGCAR format";

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);
  int VERBOSE = 0;	// The infamous verbose flag.
  int TESTING = 0;	// Extreme verbosity (testing purposes)
  int ERROR = 0;	// Analysis: Error flag (for analysis purposes)
  int MEMORY = 128;	// Our default storage (not used)
  int VORONOI = 0;	// find Voronoi volumes instead?
  int ATOM_ASSIGN = 1;	// default: assign basins to atoms
  int SCALE = 0;	// scale to the total volume of cell?
  int OUTPUT = 0;	// if -1, output all; else, index of basin to output
  char* OUTBASE = NULL;	// basename for output
  
  char ch;
  while ((ch = getopt(argc, argv, "Vaso:n:hvtm:")) != -1) {
    switch (ch) {
    case 'V':
      VORONOI = 1;
      break;
    case 'a':
      ATOM_ASSIGN = 0;
      break;
    case 's':
      SCALE = 1;
      break;
    case 'n':
      OUTPUT = (int)strtol(optarg, (char**)NULL, 10);
      break;
    case 'o':
      OUTBASE = optarg;	// should suffice...
      break;
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
  
  // if we specified an output name, but not a specific basin, output all
  if ( (OUTPUT == 0) && (OUTBASE != NULL) ) OUTPUT = -1;
  if ( (OUTPUT != 0) && (OUTBASE == NULL) ) {
    fprintf(stderr, "Need to specify output basename\n");
    ERROR = 4;
  }
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
  double** grid = new double*[argc];
  char** name = new char*[argc];
  // we have to do this duplication as basename uses static storage to do its work, and can modify the string its passed
  for (int n=0; n<argc; ++n) {
    char* dump = strdup(argv[n]);
    name[n] = strdup(basename(dump));
    delete[] dump;
  }
  if (TESTING) {
    fprintf(stderr, "## input files:\n");
    for (int n=0; n<argc; ++n)
      fprintf(stderr, "## %s\t%s\n", argv[n], name[n]);
  }

  // check input:
  FILE* infile = myopenr(infile_name);
  if (infile == NULL) {
    fprintf(stderr, "Could not open %s\n", infile_name);
    ERROR = 2;
  }

  // If we've thrown an error already, then get out.
  if (ERROR) exit(ERROR);
  if (TESTING) fprintf(stderr, "## reading %s\n", infile_name);

  // parse that file...
  char dump[512];
  char comment[512];
  double a0, alatt[9], gridlatt[9], metric[9];
  
  int* Natom=NULL;
  double** uatom;
  int Ngrid[3];

  read_CHGCAR_header(infile, comment, a0, alatt, Natom, uatom);
  square(alatt, metric);
  for (int d=0; d<9; ++d) metric[d] *= a0*a0;
  int Natoms=0;
  for (int nt=0; Natom[nt]>=0; ++nt) Natoms += Natom[nt];
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
      gridlatt[3*d+i] = a0 * alatt[3*d+i]/(double)Ngrid[i];

  if (TESTING) {
    for (int i=0; i<3; ++i) {
      printf("## A%d=", i+1);
      for (int d=0; d<3; ++d)
	printf(" %.12g", gridlatt[3*d+i]);
      printf("\n");
    }
  }

  // read the first grid:
  int Ng = Ngrid[0]*Ngrid[1]*Ngrid[2];
  double* rho = new double[Ng];
  if (!VORONOI)
    read_grid(infile, rho, Ng);
  else {
    if (TESTING) fprintf(stderr, "## generating distances for Voronoi\n");
    int iv[3], n=0;
    double uvect[3];
    for (iv[2]=0; iv[2]<Ngrid[2]; ++(iv[2])) {
      if (TESTING) printf("## %d/%d\n", iv[2], Ngrid[2]);
      uvect[2]=(double)iv[2]/(double)Ngrid[2];
      for (iv[1]=0; iv[1]<Ngrid[1]; ++(iv[1])) {
	uvect[1]=(double)iv[1]/(double)Ngrid[1];
	for (iv[0]=0; iv[0]<Ngrid[0]; ++(iv[0])) {
	  uvect[0]=(double)iv[0]/(double)Ngrid[0];
	  double r2;
	  match_atom(Natoms, uatom, metric, uvect, r2);
	  rho[n] = -r2;
	  ++n;
	}
      }
    }
    if (TESTING) fprintf(stderr, "## ... finished\n");
  }
  myclose(infile);

  if (TESTING) fprintf(stderr, "## sorting rho\n");
  // 1. Sort density in descending order
  // so that rho[order[0]] is the largest and rho[order[Ng-1]] is the smallest
  int* order = new int[Ng];
  for (int n=0; n<Ng; ++n) order[n] = n;
  // NOTE: This would be the call if we have BSD... but this version of qsort
  // doesn't appear to be general
  //  qsort_r(order, Ng, sizeof(int), rho, density_comp);
  DENSITY_COMP = rho;
  qsort(order, Ng, sizeof(int), density_comp);
  // inverse mapping, to make neighbor construction straightforward, and for
  // outputing w on a grid...
  int* redro = new int[Ng];
  for (int n=0; n<Ng; ++n) redro[order[n]]=n;
  if (TESTING) fprintf(stderr, "## ... finished\n");

  // copy in the (sorted) list
  grid[0] = new double[Ng];
#ifdef NOSORT_STORAGE
  if (!VORONOI) for (int ng=0; ng<Ng; ++ng) grid[0][ng] = rho[ng];
#else
  // or... grid[0][redro[ng]] = rho[ng]... which is how we will read the rest in.
  if (!VORONOI) for (int ng=0; ng<Ng; ++ng) grid[0][ng] = rho[order[ng]];
#endif
  else for (int ng=0; ng<Ng; ++ng) grid[0][ng] = 1;

  // now read the remaining grids...
  for (int g=1; (g<argc) && !ERROR; ++g) {      
    if (TESTING) fprintf(stderr, "## reading %s\n", argv[g]);
    infile = myopenr(argv[g]);
    if (infile == NULL) {
      fprintf(stderr, "Could not open %s\n", argv[g]);
      ERROR = 2;
      continue;
    }
    skip_CHGCAR_header(infile);
    fgets(dump, sizeof(dump), infile);
    int Ntemp[3];
    sscanf(dump, "%d %d %d", Ntemp, Ntemp+1, Ntemp+2);
    if ( (Ngrid[0] != Ntemp[0]) ||
	 (Ngrid[1] != Ntemp[1]) ||
	 (Ngrid[2] != Ntemp[2]) ) {
      fprintf(stderr, "%s does not have the same grid dimensions? %d x %d x %d\n",
	      argv[g], Ntemp[0], Ntemp[1], Ntemp[2]);
      ERROR = 4;
      continue;
    }
    grid[g] = new double[Ng];
#ifdef NOSORT_STORAGE
    read_grid(infile, grid[g], Ng);
#else
    read_grid(infile, grid[g], Ng, redro);
#endif
    myclose(infile);
  }
  if (TESTING) fprintf(stderr, "## ... finished\n");
  if (ERROR) exit(ERROR);

  // ***************************** ANALYSIS **************************
  // 2. Generate (pruned) list of neighboring vectors that bound the WS cell
  int Nneighvect;
  int** neighvect;
  double* alpha;
  
  gen_WS_voronoi(gridlatt, Nneighvect, neighvect, alpha);
  if (TESTING) {
    fprintf(stderr, "## %d\t#Nneigh\n", Nneighvect);
    for (int n=0; n<Nneighvect; ++n)
      fprintf(stderr, "## %3d%3d%3d %.15le\n",
	      neighvect[n][0], neighvect[n][1], neighvect[n][2],
	      alpha[n]);
    // fprintf(stderr, " # %.8g %.8g %.8g\n", R[n][0], R[n][1], R[n][2]);
  }

  double gridvol, volscale;
  gridvol = fabs(det(gridlatt));
  volscale = gridvol;
  if (SCALE) {
    gridvol = 1./(double)Ng;
    if (VORONOI) for (int ng=0; ng<Ng; ++ng) grid[0][ng] = volscale/gridvol;
  }
  
  // 3. Determine neighbor list and basins
  // Note: basin==0 -> interior point; basin>1 -> known to belong to a specific basin
  if (TESTING) fprintf(stderr, "## getting neighbor list and basins\n");
  int Nbasin=0;
  if (ATOM_ASSIGN) Nbasin = Natoms;
  int* basin = new int[Ng], * numbelow = new int[Ng];
  int** neigh = new int*[Ng];
  double** prob = new double*[Ng];
  // new, more efficient allocation
  neigh[0] = new int[Nneighvect*Ng];
  prob[0] = new double[Nneighvect*Ng];
  for (int ng=1; ng<Ng; ++ng) {
    neigh[ng] = neigh[ng-1] + Nneighvect;
    prob[ng] = prob[ng-1] + Nneighvect;
  }
  //  double* rho = grid[0];
  for (int n=0; n<Ng; ++n) {
    basin[n]=0;
    numbelow[n]=0;
    int i=order[n];
    int iv[3], jv[3];
    tovect(i, Ngrid, iv);
    int nabove=0;
    double rho0=rho[i];
    double t[Nneighvect], tsum=0;
    int above[Nneighvect];
    for (int nv=0; nv<Nneighvect; ++nv) {
      for (int d=0; d<3; ++d) jv[d] = iv[d] + neighvect[nv][d];
      int j=toindex(Ngrid, jv);
      int m=redro[j];
      // never though I'd have to do this, but *just in case rho[j]==rho0* ...
      //      if ( (m<n) && (rho[j]!=rho0) ) {
      if (m<n) {
	// then rho[j]>rho[i]...
	above[nabove]=m;
	t[nabove]=alpha[nv]*(rho[j]-rho0);
	tsum += t[nabove];
	++nabove;
      }
    }
    if (nabove==0) {
      // new basin!
      // we have two options... count the basins, or assign to atoms.
      if (ATOM_ASSIGN) {
	double uvect[3] = {(double)iv[0]/(double)Ngrid[0],
			   (double)iv[1]/(double)Ngrid[1],
			   (double)iv[2]/(double)Ngrid[2]};
	basin[n] = match_atom(Natoms, uatom, metric, uvect) + 1;
      } else {
	basin[n] = Nbasin + 1;
	++Nbasin;
      }
      continue;
    }
    // else, either an interior point, or a boundary point:
    //   interior == all points with larger density ("above") belong to the same basin
    //   *and* that basin != 0 (basin==0 -> boundary point)
    int tbasin = basin[above[0]];
    int boundary=0;
    for (int nneigh=0; nneigh<nabove; ++nneigh)
      if (basin[above[nneigh]] != tbasin) boundary = 1;
    boundary = boundary || (tbasin==0);
    if (boundary) {
      // we need to do some paperwork...
      basin[n]=0;
      // Error check (should very rarely be true--but when it is, it kills it all)
      if (tsum==0) {
	if (TESTING)
	  fprintf(stderr, "## found boundary point with tsum=%.6g; nabove=%d at (%d,%d,%d)\n",
	  tsum, nabove, iv[0], iv[1], iv[2]);
	for (int nneigh=0; nneigh<nabove; ++nneigh)
	  t[nneigh]=1;
	tsum=nabove;
      }
      for (int nneigh=0; nneigh<nabove; ++nneigh) {
	int m=above[nneigh];
	neigh[m][numbelow[m]] = n;
	prob[m][numbelow[m]] = t[nneigh]/tsum;
	++(numbelow[m]);
      }
    } else 
      // interior point
      basin[n] = tbasin;
  }
  if (TESTING) {
    int nboundary=0;
    for (int ng=0; ng<Ng; ++ng)
      if (basin[ng]==0) ++nboundary;
    printf("## %d / %d = %.5lf boundary\n", nboundary, Ng, 
	   (double)nboundary/(double)Ng);
  }
  if (TESTING) fprintf(stderr, "## ... finished\n");

  // 4. Loop over basins, calculating w and integrating
  printf("#BAS");
  for (int g=0; g<argc; ++g) printf(" #%s", name[g]);
  printf(" #VOL\n");

  // a little garbage collection to help keep memory reasonable
  //  delete[] rho;
  //  double* w=new double[Ng];
  double* w=rho;	// just overwrite it, rather than alloc & dealloc
  for (int nbasin=1; nbasin<=Nbasin; ++nbasin) {
    double integral[argc];
    double volume=0;
    for (int g=0; g<argc; ++g) integral[g]=0;
    for (int ng=0; ng<Ng; ++ng)
      if (basin[ng]==nbasin) w[ng]=1; else w[ng]=0;
    // now, integrate (and compute w)
    for (int ng=0; ng<Ng; ++ng) {
      double tw=w[ng];
      if (tw!=0) {
	// int* tn = neigh[ng];
	// double* tp = prob[ng];
	for (int nb=0; nb<numbelow[ng]; ++nb)
	  w[neigh[ng][nb]] += prob[ng][nb] * tw;
#ifdef NEAR_GRID
	if (tw<0.5) tw = 0; else tw = 1;
#endif
#ifdef NOSORT_STORAGE
	for (int g=0; g<argc; ++g) integral[g] += grid[g][order[ng]]*tw;
#else
	for (int g=0; g<argc; ++g) integral[g] += grid[g][ng]*tw;
#endif
	volume += tw;
      }
    }
    printf("%4d", nbasin);
    for (int g=0; g<argc; ++g) printf(" %.8le", integral[g] * gridvol);
    printf(" %.8le\n", volume*volscale);
    // code here to output weights for plotting...
    if ( (OUTPUT!=0) && ((OUTPUT == -1) || (OUTPUT == nbasin)) ) {
      sprintf(dump, "%s%04d", OUTBASE, nbasin);
      FILE* outfile = myopenw(dump);
      sprintf(dump, "%s basin %d volume %.5lf", comment, nbasin, volume);
      write_CHGCAR_header(outfile, dump, a0, alatt, Natom, uatom);
      fprintf(outfile, " %4d %4d %4d\n", Ngrid[0], Ngrid[1], Ngrid[2]);
      // note: we use redro to map from spatial index to sorted index (in w)
      write_grid(outfile, w, Ng, redro);
      myclose(outfile);
    }
  }
  delete[] w;
  
  // ************************* GARBAGE COLLECTION ********************
  delete[] basin;
  delete[] numbelow;
  delete[] neigh[0];  delete[] neigh;
  delete[] prob[0];  delete[] prob;
  
  delete[] order;
  delete[] redro;

  for (int n=0; n<Nneighvect; ++n) delete[] neighvect[n];
  delete[] neighvect;
  delete[] alpha;
  
  for (int g=0; g<argc; ++g) delete[] grid[g];
  delete[] grid;
  delete[] name;

  for (int na=0; na<Natoms; ++na) delete[] uatom[na];
  delete[] uatom;
  delete[] Natom;

  return 0;
}

// very brute force--just runs through all possible atoms looking for the shortest
// distance, and returns that atom index.    
int match_atom(const int& Natoms, double** uatom, 
	       double metric[9], double uvect[3],
	       double& closest)
{
  int best;
  double udiff[3], utry[3];
  int pbc[3];
  
  closest = metric[0] + metric[4] + metric[8];
  for (int na=0; na<Natoms; ++na) {
    for (int d=0; d<3; ++d) udiff[d] = uvect[d] - uatom[na][d];
    for (pbc[0]=-1; pbc[0]<=1; ++pbc[0]) {
      utry[0] = udiff[0] + pbc[0];
      for (pbc[1]=-1; pbc[1]<=1; ++pbc[1]) {
	utry[1] = udiff[1] + pbc[1];
	for (pbc[2]=-1; pbc[2]<=1; ++pbc[2]) {
	  utry[2] = udiff[2] + pbc[2];
	  double r2 = innerprod(utry, metric, utry);
	  if (r2 < closest) {
	    closest = r2;
	    best = na;
	  }
	}
      }
    }
  }
  return best;
}

	    
