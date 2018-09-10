/*
  Program: chgcar-transl.C
  Author:  D. Trinkle
  Date:    2006 January 10
  Purpose: We take a CHGCAR file and apply a uniform translation to the supercell.
           This involves (1) shifting atom positions, and (2) shifting the charge
	   grids.  The latter is done using FT, and does not change the grid
	   dimensions; use chgcar-mult to do that.

  Param.:  <CHGCAR> <R1> <R2> <R3>
           CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar
	   R:      displacement (cartesian, unless -u is used)

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We read in everything, and then:
           1. Multiply cell, and output, along with positions
	   2. Read the existing density grid, scale, and output
	   3. possibly repeat 2 if spin polarized

  Output:  New CHGCAR file for strained cell
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "3dfft.H"
#include "chgcar.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<CHGCAR> <R1> <R2> <R3>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'u'}; // flag characters.

const char* ARGEXPL = 
"  CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  R:      displacement (cartesian, unless -u is used)\n\
  -u      unit cell dimensions, not cartesian displacement";

const char* FILEEXPL =
"\n";

int main ( int argc, char **argv ) 
{
  int d, i; // General counting variables.

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
  int UNITCELL = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  char* chgcar_name = args[0];

  double R[3], udisp[3];
  
  sscanf(args[1], "%lf", R);
  sscanf(args[2], "%lf", R+1);
  sscanf(args[3], "%lf", R+2);
  

  infile = myopenr(chgcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_name);
    exit(ERROR_NOFILE);
  }


  // ***************************** ANALYSIS **************************
  // 0. read in what you got
  //++ comment line
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ a0
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ lattice vectors
  double inlatt[9];
  for (i=0; i<3; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", inlatt+i, inlatt+3+i, inlatt+6+i);
    printf("%s", dump);
  }

  if (! UNITCELL) {
    double invlatt[9];
    careful_inverse(inlatt, invlatt);
    mult_vect(invlatt, R, udisp);
  }
  else 
    for (d=0; d<3; ++d) udisp[d] = R[d];
  for (d=0; d<3; ++d) udisp[d] = insidecell(udisp[d]);
  
  //++ number of atoms
  fgets(dump, sizeof(dump), infile);
  char *strp, *endp=dump;
  int Natoms = 0;
  do {
    strp=endp;
    i = strtol(strp, &endp, 10);
    if (strp!=endp) Natoms += i;
  } while (strp!=endp);
  printf("%s", dump);

  //++ direct coord. (CHGCAR should always be output with "Direct")
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ atom positions
  for (i=0; i<Natoms; ++i) {
    double in_vect[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", in_vect, in_vect+1, in_vect+2);
    for (d=0; d<3; ++d) in_vect[d] = insidecell(in_vect[d] - udisp[d]);
    printf("%10.6lf%10.6lf%10.6lf\n", in_vect[0], in_vect[1], in_vect[2]);
  }
  
  //++ blank line
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);
  
  for (int ISPIN=0; ISPIN<2; ++ISPIN) {
    //++ grid values
    int N[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d %d %d", N+0, N+1, N+2);
    
    int Ng=N[0]*N[1]*N[2];
    double* rho = new double[Ng];

    read_grid(infile, rho, Ng);
    // translate!
    // 0. forward FT
    double* rhoR = new double[Ng];
    double* rhoI = new double[Ng];
    forward_3dfft(rho, N[0], N[1], N[2], rhoR, rhoI);
    // 1. translate
    int n[3], nind=0;
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
      for (n[1]=0; n[1]<N[1]; ++(n[1]))
        for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	  // phase factor: exp(iG.R')
	  double phase = 2.*M_PI*dot(n, udisp);
	  double cosP=cos(phase), sinP=sin(phase);
	  double REAL=rhoR[nind], IMAG=rhoI[nind];
	  rhoR[nind] = REAL*cosP - IMAG*sinP;
	  rhoI[nind] = REAL*sinP + IMAG*cosP;
	  ++nind;
	}

    // 2. inverse FT, and output
    inverse_3dfft(rhoR, rhoI, N[0], N[1], N[2], rho);

    // *** OUTPUT ***
    printf(" %4d %4d %4d\n", N[0], N[1], N[2]);
    write_grid(stdout, rho, Ng);
    
    // garbage collection
    delete[] rho;
    delete[] rhoR;
    delete[] rhoI;
      
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
        printf("%s", dump);
        double* aug = new double[aug_grid];
        read_grid(infile, aug, aug_grid);
        write_grid(stdout, aug, aug_grid);
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
    if (i!=0) read_grid(infile, magmom+5, i);
    write_grid(stdout, magmom, Natoms);

    delete[] magmom;
  }
  
  myclose(infile);

  return ERROR;
}
