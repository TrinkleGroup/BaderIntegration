/*
  Program: chgcar-mult.C
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: Using FFT, we'll take a CHGCAR file, read it's grid, and transform
           the file to a supercell with a new grid.

  Param.:  <CHGCAR> <super> <N1> <N2> <N3>
           CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar
	   super:  supercell definition (matrix)
	   Ni:     grid spacing for output

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We read in everything, and then:
           1. Multiply cell and positions, and output
	   2. Read the existing density grid, FFT, expand to new FFT grid, iFFT
	   3. Output augmentation charges (if they exist), multiplied
	   4. possibly repeat 2+3 if spin polarized

  Output:  New CHGCAR file for supercell on new real-space grid.
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

void get_minmax(int supern[9], int nmin[3], int nmax[3]);

// This is the heart of our algorithm, such as it is:
// returns the total density in the new cell
double conv_grid(double* rho, int N[3], 
		 double* rho_super, int M[3],
		 int super[9]);


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 5;
const char* ARGLIST = "<CHGCAR> <super> <N1> <N2> <N3>";

const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'l'}; // flag characters.

const char* ARGEXPL = 
"  CHGCAR: output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  super:  supercell definition (matrix)\n\
  Ni:     grid spacing for output\n\
  -l      explicitely allow left-handed supercells";

const char* FILEEXPL =
"==== super ====\n\
n[1,1] n[1,2] n[1,3]  n[2,1] n[2,2] n[2,3]  n[3,1] n[3,2] n[3,3]\n\
==== super ====\n";

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
  int LEFTY_OK = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  char* chgcar_name = args[0];
  char* super_name = args[1];
  
  int M[3], Mg;
  for (i=0; i<3; ++i)
    sscanf(args[i+2], "%d", &(M[i]));
  Mg = M[0]*M[1]*M[2];

  infile = myopenr(super_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", super_name);
    exit(ERROR_NOFILE);
  }
  int super[9];
  
  nextnoncomment(infile, dump, sizeof(dump));
  sscanf(dump, "%d %d %d %d %d %d %d %d %d",
	 super+0, super+1, super+2,
	 super+3, super+4, super+5,
	 super+6, super+7, super+8);
  myclose(infile);

  // invert...
  int inv_super[9];
  int Nmult = inverse(super, inv_super); 
  if (LEFTY_OK) Nmult = abs(Nmult);
  ERROR = (Nmult <= 0);

  if (ERROR) {
    fprintf(stderr, "det(super)= %d\n", Nmult);
    exit(ERROR);
  }

  infile = myopenr(chgcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_name);
    exit(ERROR_NOFILE);
  }


  // ***************************** ANALYSIS **************************
  if (TESTING) {
    printf("## Nmult= %d\n", Nmult);
    printf("## | %3d %3d %3d|\n", super[0], super[1], super[2]);
    printf("## | %3d %3d %3d|\n", super[3], super[4], super[5]);
    printf("## | %3d %3d %3d|\n", super[6], super[7], super[8]);
  }
  
  // 0. make translations for atom positions
  int Ntrans[Nmult][3];
  int n[3], nmin[3], nmax[3];
  get_minmax(super, nmin, nmax);
  if (det(super)<0) {
    for (d=0; d<3; ++d) {
      int temp = nmin[d];
      nmin[d] = -nmax[d];
      nmax[d] = -temp;
    }
  }

  if (TESTING) {
    printf("## nmin= %3d %3d %3d\n", nmin[0], nmin[1], nmin[2]);
    printf("## nmax= %3d %3d %3d\n", nmax[0], nmax[1], nmax[2]);
  }
  
  i=0;
  int superXn[3];
  for (n[0]=nmin[0]; n[0]<=nmax[0]; ++(n[0]))
    for (n[1]=nmin[1]; n[1]<=nmax[1]; ++(n[1]))
      for (n[2]=nmin[2]; n[2]<=nmax[2]; ++(n[2])) {
	mult_vect(inv_super, n, superXn);
	if ( (superXn[0] >= 0) && (superXn[0] < Nmult) &&
	     (superXn[1] >= 0) && (superXn[1] < Nmult) &&
	     (superXn[2] >= 0) && (superXn[2] < Nmult) ) {
	  for (j=0; j<3; ++j) Ntrans[i][j] = n[j];
	  ++i;
	}
      }
  if (TESTING) {
    for (j=0; j<i; ++j) 
      printf("## trans %3d = %3d %3d %3d\n", j+1,
	     Ntrans[j][0], Ntrans[j][1], Ntrans[j][2]);
    printf("##\n");
  }

  if (i != Nmult) {
    fprintf(stderr, "Somehow didn't produce the correct number of translations.\n");
    exit(1);
  }

  double unit_conv[9];
  for (d=0; d<9; ++d) unit_conv[d] = inv_super[d]/(double)Nmult;

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
  mult(inlatt, super, outlatt);
  for (i=0; i<3; ++i) 
    printf(" %12.6lf%12.6lf%12.6lf\n", outlatt[i], outlatt[i+3], outlatt[i+6]);

  //++ number of atoms
  fgets(dump, sizeof(dump), infile);
  char *strp, *endp=dump;
  int Natoms = 0;
  do {
    strp=endp;
    i = strtol(strp, &endp, 10);
    if (strp!=endp) {
      printf(" %3d", i*Nmult);
      Natoms += i;
    }
  } while (strp!=endp);
  printf("\n");

  //++ direct coord. (CHGCAR should always be output with "Direct")
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ atom positions
  double in_vect[3], temp_vect[3], new_vect[3];
  for (i=0; i<Natoms; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", in_vect, in_vect+1, in_vect+2);
    for (j=0; j<Nmult; ++j) {
      for (d=0; d<3; ++d) temp_vect[d]=in_vect[d]+Ntrans[j][d];
      mult_vect(unit_conv, temp_vect, new_vect);
      for (d=0; d<3; ++d) {
	if (new_vect[d]<0) new_vect[d] += 1.;
	if (new_vect[d]>=1) new_vect[d] -= 1.;
	if (new_vect[d]>0.999999) new_vect[d] = 0;
      }
      printf("%10.6lf%10.6lf%10.6lf\n", new_vect[0], new_vect[1], new_vect[2]);
    }
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
    double* rho_super = new double[Mg];
    
    read_grid(infile, rho, Ng);
    double total_dens = 
      conv_grid(rho, N, rho_super, M, super);
    
    printf(" %4d %4d %4d\n", M[0], M[1], M[2]);
    write_grid(stdout, rho_super, Mg);
    
    if (VERBOSE) {
      fprintf(stderr, "# average density= %.8le\n", total_dens);
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
      int nc=0; // counter for which atom we're on...
      for (int n=0; n<Natoms; ++n) {
	int aug_grid;
	sscanf(dump, "%*s %*s %*d %d", &aug_grid);
	double* aug = new double[aug_grid];
	read_grid(infile, aug, aug_grid);
	
	for (j=0; j<Nmult; ++j) {
	  ++nc;
	  printf("augmentation occupancies %3d %3d\n", nc, aug_grid);
	  write_grid(stdout, aug, aug_grid);
	}
	fgets(dump, sizeof(dump), infile);
      }
    }
    
    if (feof(infile)) {
      myclose(infile);
      return ERROR;
    }
  
    //++ we're spin-polarized, so we have to read the magmom grid, output, and go back.
    double* magmom = new double[Natoms];
    double* magmom_mult = new double[Natoms*Nmult];
    
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
    for (int n=0; n<Natoms; ++n)
      for (i=0; i<Nmult; ++i)
	magmom_mult[i+Nmult*n] = magmom[n];
    write_grid(stdout, magmom_mult, Natoms*Nmult);

    delete[] magmom;
    delete[] magmom_mult;
  }
  
  myclose(infile);

  return ERROR;
}


void get_minmax(int supern[9], int nmin[3], int nmax[3]) 
{
  int i;
  int n[3], sXn[3];
  
  for (i=0; i<3; ++i) {
    nmin[i] = 1000000;
    nmax[i] = -1000000;
  }
  
  for (n[0]=0; n[0]<2; ++(n[0]))
    for (n[1]=0; n[1]<2; ++(n[1]))
      for (n[2]=0; n[2]<2; ++(n[2])) {
	mult_vect(supern, n, sXn);
	for (i=0; i<3; ++i) {
	  if (sXn[i]<nmin[i]) nmin[i] = sXn[i];
	  if (sXn[i]>nmax[i]) nmax[i] = sXn[i];
	}
      }
}



// This is the heart of our algorithm, such as it is:
double conv_grid(double* rho, int N[3], 
		 double* rho_super, int M[3],
		 int super[9]) 
{
  // 0. FFT
  long Ng = N[0]*N[1]*N[2];
  double *rhoR = new double[Ng], *rhoI = new double[Ng];
  
  forward_3dfft(rho, N[0], N[1], N[2], rhoR, rhoI);
  
  // 1. setup new grid
  long Mg = M[0]*M[1]*M[2];
  double *rhoSR = new double[Mg], *rhoSI = new double[Mg];
  for (long n=0; n<Mg; ++n) {
    rhoSR[n]=0;
    rhoSI[n]=0;
  }

  // 2. copy equivalent G's in original rlv to new G's in new rlv
  double scale = det(super)*Mg / (double)Ng;
  if (scale < 0) scale = -scale;
  int supert[9];
  transpose(super, supert);
  if (det(super) < 0) mult(supert, -1, supert);
  long M0 = M[0], M01=M[0]*M[1];

  int n[3], nc[3], nc_h[3], nnew[3];
  int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
  
  long nind=0, nind_new;
  int h[3], half[3];
  for (n[2]=0; n[2]<N[2]; ++(n[2]))
    for (n[1]=0; n[1]<N[1]; ++(n[1]))
      for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	for (int d=0; d<3; ++d) {
	  if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
	  else            nc[d] = n[d];
	  if (n[d]==Nh[d]) half[d] = 1;
	  else             half[d] = 0;
	}
	double hscale = scale/( (half[0]+1)*(half[1]+1)*(half[2]+1) );
	double scaled_rhoR = hscale*rhoR[nind];
	double scaled_rhoI = hscale*rhoI[nind];
	++nind;
	
	for (h[0]=0; h[0]<=half[0]; ++(h[0]))
	  for (h[1]=0; h[1]<=half[1]; ++(h[1]))
	    for (h[2]=0; h[2]<=half[2]; ++(h[2])) {
	      for (int d=0; d<3; ++d) nc_h[d] = nc[d] - h[d]*N[d];
	      mult_vect(supert, nc_h, nnew);
	      for (int d=0; d<3; ++d) nnew[d] = (nnew[d]%M[d]+M[d])%M[d];
	      nind_new = nnew[0] + nnew[1]*M0 + nnew[2]*M01;
	      rhoSR[nind_new] += scaled_rhoR;
	      rhoSI[nind_new] += scaled_rhoI;
	    }
      }

  double total_dens = rhoSR[0]/Mg;
  
  // 3. iFFT
  int not_real =
    inverse_3dfft(rhoSR, rhoSI, M[0], M[1], M[2], rho_super);
  if (not_real) {
    fprintf(stderr, "Somehow... did not produce a real density?\n");
  }

  // 4. Garbage collection
  delete[] rhoR;  delete[] rhoI;
  delete[] rhoSR; delete[] rhoSI;

  return total_dens;
}
