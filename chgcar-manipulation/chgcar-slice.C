/*
  Program: chgcar-slice.C
  Author:  D. Trinkle
  Date:    2010 December 3
  Purpose: Read a CHGCAR-formated file and output a planar slice.

           Right now, this is VERY simple--it has to be a grid slice.  However,
	   with some modification it should be possible to do a general
	   planar slice, using tricubic spline package.

  Param.:  CHGCAR plane u_plane
           CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar
	   plane:  integer index (1,2,3) for which slice to take
	   u:      displacement (between 0 and 1)

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   

  Output:  Grid of xy positions and values.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <libgen.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "3dfft.H"
#include "chgcar.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 3;
const char* ARGLIST = "[-s] [-m MEMORY] [-hvt] CHGCAR plane u_plane";

const char* ARGEXPL = 
"  CHGCAR:	output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  plane:	integer index (1,2,3) for which slice to take\n\
  u:	displacement (between 0 and 1)\n\
  -s	scale by volume to output a density (charge/A^3)\n\
\n\
  -m [MEMORY]  suggest necessary memory for input\n\
  -h  HELP\n\
  -v  VERBOSE\n\
  -t  TESTING";

const char* FILEEXPL =
"VASP CHGCAR (real-space grid) format\n";

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter
  int SCALE = 0;

  char ch;
  while ((ch = getopt(argc, argv, "shvtm:")) != -1) {
    switch (ch) {
    case 's':
      SCALE = -1;
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

  if (TESTING) {
    printf("# MEMORY=%d\n", MEMORY);
  }
  if (VERBOSE) {
  }

  if (argc < 1) ERROR = 2;
  // All hell broken loose yet?
  if (ERROR != 0) {
    fprintf(stderr, "%s %s\n", progname, ARGLIST);
    fprintf(stderr, "%s\n", ARGEXPL);
    fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
    exit(ERROR);
  }

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  char* chgcar_name = argv[0];
  int plane = strtol(argv[1], NULL, 10);
  double u_plane = strtod(argv[2], NULL);
  plane--;	// so that 1->0, 2->1, 3->2 for indexing
  u_plane -= (int)(u_plane);
  if (u_plane<0) u_plane += 1.;

  infile = myopenr(chgcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_name);
    exit(ERROR_NOFILE);
  }

  //++ comment line
  fgets(dump, sizeof(dump), infile);

  //++ a0
  double a0;
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%lf", &a0);

  //++ lattice vectors
  double inlatt[9], outlatt[9];
  for (int i=0; i<3; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", inlatt+i, inlatt+3+i, inlatt+6+i);
  }
  for (int i=0; i<9; ++i) inlatt[i] *= a0;
  square(inlatt, outlatt);
  double vol = det(inlatt);
  double len[3] = {sqrt(outlatt[0]), sqrt(outlatt[4]), sqrt(outlatt[8])};

  //++ number of atoms
  fgets(dump, sizeof(dump), infile);
  char *strp, *endp=dump;
  int Natoms = 0;
  do {
    strp=endp;
    int i = strtol(strp, &endp, 10);
    if (strp!=endp) Natoms += i;
  } while (strp!=endp);

  //++ direct coord. (CHGCAR should always be output with "Direct")
  fgets(dump, sizeof(dump), infile);

  //++ atom positions
  double in_vect[3];
  for (int i=0; i<Natoms; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", in_vect, in_vect+1, in_vect+2);
  }
  
  //++ blank line
  fgets(dump, sizeof(dump), infile);
  
  // ***************************** ANALYSIS **************************

  for (int ISPIN=0; ISPIN<2; ++ISPIN) {
    //++ grid values
    int N[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d %d %d", N+0, N+1, N+2);
    
    int Ng=N[0]*N[1]*N[2];
    double* rho = new double[Ng];
    
    read_grid(infile, rho, Ng);
    
    int u[3];
    int u_match = ((int)(u_plane*N[plane]))%(N[plane]);
    int ind=0;
    int planex = (plane+1)%3, planey = (plane+2)%3;
    double dX = len[planex]/N[planex];
    double dY = len[planey]/N[planey];
    double scale = 1;
    if (SCALE) scale = 1./vol;
    for (u[2]=0; u[2]<N[2]; ++u[2]) 
      for (u[1]=0; u[1]<N[1]; ++u[1])
	for (u[0]=0; u[0]<N[0]; ++u[0]) {
	  if (u[plane] == u_match) {
	    double x = u[planex]*dX;
	    double y = u[planey]*dY;
	    printf("%.6lf %.6lf %.12g\n", x, y, rho[ind]*scale);
	  }
	  ++ind;
	}
    
    // garbage collection
    delete[] rho;
    
    if (feof(infile)) {
      myclose(infile);
      return ERROR;
    }
    
    //++ (augmentation charges?)
    fgets(dump, sizeof(dump), infile);
    if (dump[0] == 'a') {
      // augmentation charges...
      for (int n=0; n<Natoms; ++n) {
	int aug_grid;
	sscanf(dump, "%*s %*s %*d %d", &aug_grid);
	double* aug = new double[aug_grid];
	read_grid(infile, aug, aug_grid);
	fgets(dump, sizeof(dump), infile);
      }
    }
    
    if (feof(infile)) {
      myclose(infile);
      return ERROR;
    }
  
    //++ we're spin-polarized, so we have to read the magmom grid, output, and go back.
    double* magmom = new double[Natoms];
    
    // ... parse the line we've already read.
    int i=Natoms-5;
    if (i<0) i=0;
    {
      char *strp, *endp=dump;
      for (int n=0; n<(Natoms-i); ++n) {
	strp=endp;
	magmom[n] = strtod(strp, &endp);
      }
    }

    delete[] magmom;
  }
  
  myclose(infile);

  return ERROR;
}
