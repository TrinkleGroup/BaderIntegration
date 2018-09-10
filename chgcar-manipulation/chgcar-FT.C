/*
  Program: chgcar-FT.C
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: Using FFT, we'll take a CHGCAR file, read it's grid, and FT
           the file, outputting magnitude of rho(G) with |G|.

  Param.:  <CHGCAR>
           CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We read in everything, and then FT, and output. (not sorted)

  Output:  
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

void read_grid(FILE* infile, double* rho, int Ng);

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 1;
const char* ARGLIST = "<CHGCAR>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'p'}; // flag characters.

const char* ARGEXPL = 
"  CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  -p      ignore magnetization density (if present)";

const char* FILEEXPL =
"\n";

int main ( int argc, char **argv ) 
{
  int d, i, j; // General counting variables.

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
  int NOSPIN = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  char* chgcar_name = args[0];
  
  infile = myopenr(chgcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_name);
    exit(ERROR_NOFILE);
  }

  // ***************************** ANALYSIS **************************
  //++ comment line
  fgets(dump, sizeof(dump), infile);

  //++ a0
  fgets(dump, sizeof(dump), infile);

  //++ lattice vectors
  double latt[9];
  for (i=0; i<3; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", latt+i, latt+3+i, latt+6+i);
  }

  //++ number of atoms
  fgets(dump, sizeof(dump), infile);
  char *strp, *endp=dump;
  int Natoms = 0;
  do {
    strp=endp;
    i = strtol(strp, &endp, 10);
    if (strp!=endp) Natoms += i;
  } while (strp!=endp);

  //++ direct coord. (CHGCAR should always be output with "Direct")
  fgets(dump, sizeof(dump), infile);

  //++ atom positions
  for (i=0; i<Natoms; ++i) {
    fgets(dump, sizeof(dump), infile);
  }
  
  //++ blank line
  fgets(dump, sizeof(dump), infile);

  double Gmetric[9];
  {
    double temp[9];
    careful_inverse(latt, temp);	// temp = [a]^-1
    mult(temp, 2*M_PI, temp);		// temp = 2pi[a]^-1
    square(temp, Gmetric);		// Gmetric = (2pi)^2 [a]^-T[a]^-1
  }

  int MAXSPIN=2;
  if (NOSPIN) MAXSPIN=1;

  for (int ISPIN=0; ISPIN<MAXSPIN; ++ISPIN) {
    //++ grid values
    int N[3];
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%d %d %d", N+0, N+1, N+2);
    
    int Ng=N[0]*N[1]*N[2];
    double* rho = new double[Ng];
    double* rhoR = new double[Ng];
    double* rhoI = new double[Ng];
    
    read_grid(infile, rho, Ng);
    forward_3dfft(rho, N[0], N[1], N[2], rhoR, rhoI);
   
    // ************* OUTPUT ***************
    int n[3], nc[3], nc_h[3];
    int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
  
    int nind=0;
    int h[3], half[3];
    double Gabs;
    double peratom_scale = 1./((double)Natoms * (double)Ng);
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
      for (n[1]=0; n[1]<N[1]; ++(n[1]))
	for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	  for (int d=0; d<3; ++d) {
	    if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
	    else            nc[d] = n[d];
	    if (n[d]==Nh[d]) half[d] = 1;
	    else             half[d] = 0;
	  }
	  double scaled_rho = hypot(rhoR[nind], rhoI[nind]) * peratom_scale
	    /( (half[0]+1)*(half[1]+1)*(half[2]+1) );
	  ++nind;
        
	  for (h[0]=0; h[0]<=half[0]; ++(h[0]))
	    for (h[1]=0; h[1]<=half[1]; ++(h[1]))
	      for (h[2]=0; h[2]<=half[2]; ++(h[2])) {
		for (int d=0; d<3; ++d) nc_h[d] = nc[d] - h[d]*N[d];
		Gabs = sqrt(magnsq(Gmetric, nc_h));
		printf("%10.5lf %.5le\n", Gabs, scaled_rho);
	      }
	}
    
 
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
    if (i!=0) read_grid(infile, magmom+5, i);

    delete[] magmom;

    if (!NOSPIN) printf("&\n");
  }
  
  myclose(infile);

  return ERROR;
}


// this is how many entries appear per line...
const int READ_stride = 5;

void read_grid(FILE* infile, double* rho, int Ng) 
{
  int n, line, nl;
  char dump[512];

  nl = Ng/READ_stride;
  n=0;
  
  for (line=0, n=0; line<nl; ++line, n+=READ_stride) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf %lf %lf",
	   rho+n, rho+n+1, rho+n+2, rho+n+3, rho+n+4);
  }

  if (n!=Ng) {
    // last line...
    fgets(dump, sizeof(dump), infile);
    char *strp, *endp=dump;
    for ( ; n<Ng; ++n) {
      strp=endp;
      rho[n] = strtod(strp, &endp);
    }
  }
}
