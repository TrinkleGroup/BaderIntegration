/*
  Program: test-FFT.C
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: Test our implementation of the 3d FFT / iFFT code.

  Param.:  <...>

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We read in our positions and forces, and convert the position
           to unit cell coord.  The forces become dynamical matrix elements
	   after inverting the displacements.

  Output:  Usual lattice function output
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "3dfft.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 0;
const char* ARGLIST = "<...>";


const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"";

const char* FILEEXPL =
"\n";

int main ( int argc, char **argv ) 
{
  int i, j, k; // General counting variables.

  // ************************** INITIALIZATION ***********************
  int VERBOSE = 0;  // The infamous verbose flag.
  int TESTING = 0;  // Extreme verbosity (testing purposes)
  int ERROR = 0;    // Analysis: Error flag (for analysis purposes)
  int MEMORY = 1024;// Our default gridding parameter

  char* args[NUMARGS];
  int flagon[NFLAGS]; // We use this to determine which flags are on.
  for (i=0; i<NFLAGS; ++i) flagon[i] = 0; // Set all flags off.

  // Read our commandline.
  ERROR = parse_commandline(argc, argv, NUMARGS, args,
			    VERBOSE, TESTING, MEMORY, 
			    NFLAGS, USERFLAGLIST, flagon);
  // All hell broken loose yet?
  if (ERROR != 0) {
    if ( has_error(ERROR, ERROR_HELP) ) {
      // Note: we don't print out the elastic const. info... mainly
      // because we just ignore all that stuff anyway.
      print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
      fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
      fprintf(stderr, "\n");
    }
    else print_short_help(argv[0], ARGLIST, NFLAGS, USERFLAGLIST, ARGEXPL);
    exit(ERROR);
  }

  // flags

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;


  double *Z, *Zt, *Xr, *Xi;

  if (ERROR) exit(ERROR);

  // ***************************** ANALYSIS **************************

  //  int N[3] = {4,4,4}, Ng;
  int N[3] = {4,7,5}, Ng;
  int n[3];
  
  Ng=N[0]*N[1]*N[2];
  Z = new double[Ng];
  Zt = new double[Ng];
  Xr = new double[Ng];
  Xi = new double[Ng];
  
  int nind=0;
  double dN[3] = {1./(double)N[0], 1./(double)N[1], 1./(double)N[2]};
  
  for (n[2]=0; n[2]<N[2]; ++(n[2]))
    for (n[1]=0; n[1]<N[1]; ++(n[1]))
      for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	double u = n[0]*dN[0]+n[1]*dN[1]+n[2]*dN[2];
	Z[nind] = cos(2*M_PI*u);
	++nind;
      }


  forward_3dfft(Z, N[0], N[1], N[2], Xr, Xi);
  int not_real =
    inverse_3dfft(Xr, Xi, N[0], N[1], N[2], Zt);

  // ****************************** OUTPUT ***************************


  if (VERBOSE) {
    printf("Ng= %d %d %d (%d)\n", N[0], N[1], N[2], Ng);
    
    nind=0;
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
      for (n[1]=0; n[1]<N[1]; ++(n[1]))
	for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	  printf("Z[ %2d %2d %2d ]= %.8lf\n", n[0], n[1], n[2], Z[nind]);
	  ++nind;
	}
    printf("\n");
  }

  printf("Ng= %d %d %d (%d)\n", N[0], N[1], N[2], Ng);
  
  nind=0;
  for (n[2]=0; n[2]<N[2]; ++(n[2]))
    for (n[1]=0; n[1]<N[1]; ++(n[1]))
      for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	printf("X[ %2d %2d %2d ]= %.8lf + i %.8lf\n", n[0], n[1], n[2], 
	       Xr[nind], Xi[nind]);
	++nind;
      }
  printf("\n");


  if (VERBOSE) {
    printf("Ng= %d %d %d (%d)\n", N[0], N[1], N[2], Ng);
    
    nind=0;
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
      for (n[1]=0; n[1]<N[1]; ++(n[1]))
	for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	  printf("Z[ %2d %2d %2d ]= %.8lf\n", n[0], n[1], n[2], Z[nind]);
	  ++nind;
	}
    printf("\n");
  }

  printf("inverse transform back to Z: not_real= %d\n", not_real);
  double diff_sum = 0, diff;
  for (i=0; i<Ng; ++i) {
    diff = Z[i]-Zt[i];
    diff_sum += diff*diff;
  }
  diff_sum *= 1./(double)Ng;
  printf("average difference= %.8le\n", diff_sum);

  // ************************* GARBAGE COLLECTION ********************

  delete[] Z;
  delete[] Zt;
  delete[] Xr;
  delete[] Xi;

  return ERROR;
}
