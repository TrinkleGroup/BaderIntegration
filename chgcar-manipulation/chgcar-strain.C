/*
  Program: chgcar-strain.C
  Author:  D. Trinkle
  Date:    2005 December 29
  Purpose: We take a CHGCAR file and apply a strain to the supercell.  This
           involves (1) transforming the lattice, and (2) scaling the charge
	   by the total volume.  It does not change the grid, etc.  If that
	   is your goal, use chgcar-mult.

  Param.:  <CHGCAR> <strain>
           CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar
	   strain: matrix (not assumed symmetric)

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
const int NUMARGS = 2;
const char* ARGLIST = "<CHGCAR> <strain>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  strain: matrix (not assumed symmetric)";

const char* FILEEXPL =
"==== strain ====\n\
eps[x,x] eps[x,y] eps[x,z]  eps[y,x] eps[y,y] eps[y,z]  eps[z,x] eps[z,y] eps[z,z]\n\
==== strain ====\n";

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  char* chgcar_name = args[0];
  char* strain_name = args[1];
  
  infile = myopenr(strain_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", strain_name);
    exit(ERROR_NOFILE);
  }
  double strain[9];
  
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	 strain+0, strain+1, strain+2,
	 strain+3, strain+4, strain+5,
	 strain+6, strain+7, strain+8);
  myclose(infile);

  // add the identity:
  for (d=0; d<9; d += 4) strain[d] += 1;

  // invert...
  double scale = det(strain);
  ERROR = zero(scale);

  if (ERROR) {
    fprintf(stderr, "det(strain)= %.8le\n", scale);
    exit(ERROR);
  }

  infile = myopenr(chgcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_name);
    exit(ERROR_NOFILE);
  }


  // ***************************** ANALYSIS **************************
  if (TESTING) {
    printf("## scale= %.12le\n", scale);
    printf("## 1 + eps =\n");
    printf("## | %8.5lf %8.5lf %8.5lf|\n", strain[0], strain[1], strain[2]);
    printf("## | %8.5lf %8.5lf %8.5lf|\n", strain[3], strain[4], strain[5]);
    printf("## | %8.5lf %8.5lf %8.5lf|\n", strain[6], strain[7], strain[8]);
    printf("##\n");
  }

  //++ comment line
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ a0
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ lattice vectors
  double inlatt[9], outlatt[9];
  for (i=0; i<3; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", inlatt+i, inlatt+3+i, inlatt+6+i);
  }
  // apply strain!! -- this is it:
  mult(strain, inlatt, outlatt);
  for (i=0; i<3; ++i) 
    printf(" %12.6lf%12.6lf%12.6lf\n", outlatt[i], outlatt[i+3], outlatt[i+6]);

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
    fgets(dump, sizeof(dump), infile);
    printf("%s", dump);
  }
  
  //++ blank line
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);
  
  for (int ISPIN=0; ISPIN<2; ++ISPIN) {
    //++ grid values
    int N[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d %d %d", N+0, N+1, N+2);
    printf("%s", dump);
    
    int Ng=N[0]*N[1]*N[2];
    double* rho = new double[Ng];

    read_grid(infile, rho, Ng);
    // scale!
    for (int n=0; n<Ng; ++n) rho[n] *= scale;
    
    write_grid(stdout, rho, Ng);
    
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
