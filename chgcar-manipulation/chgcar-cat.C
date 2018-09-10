/*
  Program: chgcar-cat.C
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: Using linear interpolation, we'll take two CHGCAR files, read their
           grids, and then "concatenate" the two together to produce a new
	   CHGCAR.  We also smooth everything using a low-pass filter, similar
	   to what we did for the dislocation case.

  Param.:  <CHGCAR.1> <CHGCAR.2> <axis> <N1> <N2> <N3>
           CHGCAR.i:  output from VASP; handles US, PAW, spin- and non-spin-polar
	   axis:      1,2,3 -- corresponding to which axis to concatenate
	   N.i:       grid for new CHGCAR

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We read and compute as needed:           
           1. Read headers for both, and check for compatibility
	   2. Output atom positions in combined cell
	   3. Read each density grid, and compute combined low-pass filter
	   4. Construct new density by linear interpolation in each region
	      of space
	   5. FFT, apply a low-pass filter, and iFFT back; then output
	   6. Output augmentation charges (if they exist) akin to #2
	   7. repeat if spin-polarized--need to handle case where only one
	      CHGCAR is polarized by making the other have 0 spin.

  Output:  New CHGCAR file for dislocation.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "3dfft.H"
#include "linear.H"
#include "chgcar.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

void extract_filter(double* rhoR, double* rhoI, 
		    int N[3], double latt[9],
		    double& filter_a, double& filter_b);

void apply_filter(double* rhoR, double* rhoI, 
		  int N[3], double latt[9],
		  double filter_a, double filter_b);

inline int sum_atoms(int* Natom) 
{
  int sum=0;
  for (int n=0; Natom[n]>=0; ++n)
    sum+=Natom[n];
  return sum;
}

void dealloc_atoms(int Natoms, double** &uatom) 
{
  for (int n=0; n<Natoms; ++n) delete[] uatom[n];
  delete[] uatom;
  uatom=NULL;
}

int lattice_match(double latt0[9], double latt1[9], int axis, double &alpha);


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 6;
const char* ARGLIST = "<CHGCAR.1> <CHGCAR.2> <axis> <N1> <N2> <N3>";


const int NFLAGS = 1;
const char USERFLAGLIST[NFLAGS] = {'f'}; // flag characters.

const char* ARGEXPL = 
"  CHGCAR.i:  output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  axis:      1,2,3 -- corresponding to which axis to concatenate\n\
  N.i:       grid for new CHGCAR\n\
  -f         NO filtering (just concatenate--intended for testing purposes)\n\
\n\
*** WARNING: CHGCAR.1 must not be the same file as CHGCAR.2";

const char* FILEEXPL =
"NOTE: CHGCAR's MUST be compatible--\n\
  1. Equal number of atom *types*--this may necessitate adding 0 values\n\
  2. Cells that can be concatenated along the given axis.\n";

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
  int NOFILTER = flagon[0];

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile[2];

  char* chgcar1_name = args[0];
  char* chgcar2_name = args[1];
  int axis;
  int N[3], Ng;
  sscanf(args[2], "%d", &axis);
  sscanf(args[3], "%d", &(N[0]));
  sscanf(args[4], "%d", &(N[1]));
  sscanf(args[5], "%d", &(N[2]));
  Ng = N[0]*N[1]*N[2];
  if ( (axis < 1) || (3 < axis) ) {
    fprintf(stderr, "axis (%d) must be 1,2,3\n", axis);
    exit(ERROR_BADFLAG);
  }
  axis--; // set to 0,1,2
  
  if ( (N[0]<1) || (N[1]<1) || (N[2]<1) ) {
    fprintf(stderr, "grid (%d x %d x %d) must be positive\n", N[0],N[1],N[2]);
    exit(ERROR_BADFLAG);
  }
  
  infile[0] = myopenr(chgcar1_name);
  infile[1] = myopenr(chgcar2_name);
  if (infile[0] == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar1_name);
    exit(ERROR_NOFILE);
  }
  if (infile[1] == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar2_name);
    exit(ERROR_NOFILE);
  }

  // ***************************** ANALYSIS **************************
  char comment[2][512];
  double a0[2];
  double latt[2][9];
  int Natom[2][128];
  //  double*** uatom = new (double**)[2];
  double** uatom[2];
  
  for (int f=0; f<2; ++f)
    read_CHGCAR_header(infile[f], comment[f], a0[f], latt[f], Natom[f], uatom[f]);

  if (TESTING) {
    for (int f=0; f<2; ++f) {
      printf("==== CHGCAR.%d header ====\n", f+1);
      write_CHGCAR_header(stdout, comment[f], a0[f], latt[f], Natom[f], uatom[f]);
    }
  }

  // 0. Test compatibility, and output header.
  // 0.a. atom types equal
  for (i=0; (Natom[0][i]>=0) && (Natom[1][i]>=0); ++i) ;
  if ( (Natom[0][i]>=0) || (Natom[1][i]>=0) ) {
    fprintf(stderr, "Incompatible atom types?\n");
    for (int f=0; f<2; ++f) {
      fprintf(stderr, "CHGCAR.%d: ", f+1);
      for (i=0; (Natom[f][i]>=0); ++i) fprintf(stderr, " %d", Natom[f][i]);
      fprintf(stderr, "\n");
    }
    exit(ERROR_BADFILE);
  }

  // 0.b lattice constants equal
  if (! dcomp(a0[0], a0[1])) {
    fprintf(stderr, "Sorry--we can\'t handle differently scaled lattices yet.\n");
    fprintf(stderr, "a0.1= %.16lf  a0.2= %.16lf\n", a0[0], a0[1]);
    exit(ERROR_BADFILE);
  }
  
  // 0.c lattices match up
  double alpha=0;
  if (! lattice_match(latt[0], latt[1], axis, alpha)) {
    fprintf(stderr, "Sorry--we need lattice faces to match up, with parallel axes.\n");
    exit(ERROR_BADFILE);
  }

  if (VERBOSE) {
    fprintf(stderr, "# alpha= %.9lf\n", alpha);
  }

  // 0.d header output
  double lattc[9];
  int Natc[128], Nc = sum_atoms(Natom[0])+sum_atoms(Natom[1]);
  {
    char commentc[1024];
    strcpy(commentc, comment[0]);
    strcat(commentc, comment[1]);
    double** uatc = new double*[Nc];
    for (int n=0; n<Nc; ++n) uatc[n] = new double[3];
    double scale_mat[9];
    for (d=0; d<9; ++d) scale_mat[d] = ident[d];
    scale_mat[axis*4] = 1/alpha;
    mult(latt[0], scale_mat, lattc);

    int nc[3]={0,0,0}; // [0]: counter in 0, [1]: in 1, [2]: in combined
    double scale[2]={alpha, 1-alpha};
    double shift[2]={0, alpha};

    int t;
    for (t=0; Natom[0][t]>=0; ++t) {
      for (int f=0; f<2; ++f) 
	for (i=0; i<Natom[f][t]; ++i) {
	  for (d=0; d<3; ++d)
	    uatc[nc[2]][d] = uatom[f][nc[f]][d];
	  uatc[nc[2]][axis] = scale[f]*uatc[nc[2]][axis] + shift[f];
	  ++(nc[2]);
	  ++(nc[f]);
	}
      Natc[t] = Natom[0][t]+Natom[1][t];
    }
    Natc[t] = -1;
    write_CHGCAR_header(stdout, commentc, a0[0], lattc, Natc, uatc);
    dealloc_atoms(Nc, uatc);
  }

  double* rhoc;
  double scale_ind[2][3];
  
  // for (int ISPIN=0; ISPIN<1; ++ISPIN) {
  for (int ISPIN=0; ISPIN<2; ++ISPIN) {
    //++ grid values
    int M[2][3], Mg[2];
    double* rho[2];
    for (int f=0; f<2; ++f) {
      fgets(dump, sizeof(dump), infile[f]);
      sscanf(dump, "%d %d %d", M[f]+0, M[f]+1, M[f]+2);
    
      Mg[f]=M[f][0]*M[f][1]*M[f][2];
      rho[f] = new double[Mg[f]];
    
      read_grid(infile[f], rho[f], Mg[f]);
    }

    for (int f=0; f<2; ++f)
      for (d=0; d<3; ++d) scale_ind[f][d] = (double)(M[f][d])/(double)(N[d]);
    scale_ind[0][axis] *= 1/alpha;
    scale_ind[1][axis] *= 1/(1-alpha);
  
    rhoc = new double[Ng];
    
    // 0. Linear interpolation
    int n[3], nind=0;
    double u[3];
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
      for (n[1]=0; n[1]<N[1]; ++(n[1]))
	for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	  for (d=0; d<3; ++d)
	    u[d] = n[d] * scale_ind[0][d];
	  // here's where we decide if we're in region I or II
	  if (u[axis]<=M[0][axis]) 
	    rhoc[nind] = linear(rho[0], M[0], u) * (1/alpha);
	  else {
	    for (d=0; d<3; ++d)
	      u[d] = n[d] * scale_ind[1][d];
	    u[axis] = (n[axis]-alpha*N[axis])*scale_ind[1][axis];
	    rhoc[nind] = linear(rho[1], M[1], u) * (1/(1-alpha));
	  }
	  ++nind;
	}
    
    if (! NOFILTER) {
      // Now, all of our analysis.
      // 0. determine our low-pass filter values
      double filter_a[2], filter_b[2];
      for (int f=0; f<2; ++f) {
	double* rhoR = new double[Mg[f]];
	double* rhoI = new double[Mg[f]];
	forward_3dfft(rho[f], M[f][0], M[f][1], M[f][2], rhoR, rhoI);
	extract_filter(rhoR, rhoI, M[f], latt[f], filter_a[f], filter_b[f]);
	if (VERBOSE) {
	  fprintf(stderr, "# rho_filter_%d(G)= %.8le |G|^ %.8lf\n",
		  f+1, filter_a[f], filter_b[f]);
	}
	// garbage collection
	delete[] rhoR;
	delete[] rhoI;
      }

      // 1. compute combined filter value:
      const double FILTER_MULT = 3.0;  // at some point, should make this a param...
      double filtc_a = FILTER_MULT*(filter_a[0]/alpha + filter_a[1]/(1-alpha));
      double filtc_b = filter_b[0];
      if (filter_b[1]>filter_b[0]) filtc_b = filter_b[1];
      if (VERBOSE) {
	fprintf(stderr, "# rho_filt_comb(G)= %.8le |G|^ %.8lf\n", filtc_a, filtc_b);
      }
      
      // 2. low-pass truncation filter
      double* rhoR = new double[Ng];
      double* rhoI = new double[Ng];
      forward_3dfft(rhoc, N[0], N[1], N[2], rhoR, rhoI);
      apply_filter(rhoR, rhoI, N, lattc, filtc_a, filtc_b);
      inverse_3dfft(rhoR, rhoI, N[0], N[1], N[2], rhoc);
      delete[] rhoR;
      delete[] rhoI;
    }
    

    for (int f=0; f<2; ++f) delete[] rho[f];

    // *** OUTPUT ***
    printf(" %4d %4d %4d\n", N[0], N[1], N[2]);
    write_grid(stdout, rhoc, Ng);

    delete[] rhoc;

    //++ (augmentation charges?)
    //NOTE: augmentation charge lines appear to be there even for US pseudopotentials
    //... so I'm writing this ASSUMING those lines will be there.
    // The loop structure is taken from the header output above to put the atom
    // positions.
    {
      int nc[3]={0,0,0}; // our combined counter as before.
      for (int t=0; Natom[0][t]>=0; ++t) {
	for (int f=0; f<2; ++f) 
	  for (i=0; i<Natom[f][t]; ++i) {	
	    fgets(dump, sizeof(dump), infile[f]);
	    int aug_grid;
	    sscanf(dump, "%*s %*s %*d %d", &aug_grid);
	    printf("augmentation occupancies %3d %3d\n", nc[2]+1, aug_grid);
	    double* aug = new double[aug_grid];
	    read_grid(infile[f], aug, aug_grid);
	    write_grid(stdout, aug, aug_grid);
	    ++(nc[f]);
	    ++(nc[2]);
	    delete[] aug;
	  }
      }
    }
    

    if (feof(infile[0]) || feof(infile[1])) {
      for (int f=0; f<2; ++f) {
	myclose(infile[0]);
	dealloc_atoms(sum_atoms(Natom[f]), uatom[f]);
      }
      return ERROR;
    }

    //++ we're spin-polarized, so we have to read the magmom grid, output, and go back.
    double* magmom[3] = {new double[sum_atoms(Natom[0])],
			 new double[sum_atoms(Natom[1])],
			 new double[Nc]};

    for (int f=0; f<2; ++f)
      read_grid(infile[f], magmom[f], sum_atoms(Natom[f]));

    if (feof(infile[0]) || feof(infile[1])) {
      for (int f=0; f<2; ++f) {
	myclose(infile[0]);
	dealloc_atoms(sum_atoms(Natom[f]), uatom[f]);
	delete[] magmom[f];
      }
      delete[] magmom[2];
      return ERROR;
    }

    // Same loop structure again.
    {
      int nc[3]={0,0,0}; // our combined counter as before.
      for (int t=0; Natom[0][t]>=0; ++t) {
	for (int f=0; f<2; ++f) 
	  for (i=0; i<Natom[f][t]; ++i) {	
	    magmom[2][nc[2]] = magmom[f][nc[f]];
	    ++(nc[f]);
	    ++(nc[2]);
	  }
      }
      write_grid(stdout, magmom[2], Nc);
    }
    // Garbage collection
    for (int f=0; f<3; ++f) delete[] magmom[f];

    // ...and go back to the beginning to deal with spin.
  }

  for (int f=0; f<2; ++f) {
    myclose(infile[f]);
    dealloc_atoms(sum_atoms(Natom[f]), uatom[f]);
  }

  return ERROR;
}


void add_list (double* G_l, double* rho_l, int& Nl, double G, double rho) 
{
  if (zero(rho)) return;

  if (Nl==0) {
    Nl++;
    G_l[0] = G;
    rho_l[0] = rho;
    return;
  }

  int i;
  for (i=0; (i<Nl) && (G>G_l[i]); ++i) ;
  if (dcomp(G_l[i], G)) {
    if (rho_l[i] < rho) rho_l[i] = rho;
  } else {
    for (int j=Nl; j>i; --j) {
      G_l[j] = G_l[j-1];
      rho_l[j] = rho_l[j-1];
    }
    G_l[i] = G;
    rho_l[i] = rho;
    ++Nl;
  }
}


void extract_filter(double* rhoR, double* rhoI, 
		    int N[3], double latt[9],
		    double& filter_a, double& filter_b) 
{
  double Gmetric[9];
  {
    double temp[9];
    careful_inverse(latt, temp);	// temp = [a]^-1
    mult(temp, 2*M_PI, temp);		// temp = 2pi[a]^-1
    square(temp, Gmetric);		// Gmetric = (2pi)^2 [a]^-T[a]^-1
  }
    
  int n[3], nc[3], nc_h[3];
  int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
  
  int nind=0;
  int h[3], half[3];
  double Gabs;
  int Nl = 131072; // hardcoded, I know...
  double* G_l = new double[Nl];
  double* rho_l = new double[Nl];
  Nl = 0;

  double rho_min = rhoR[0] * 1e-4;
  
  for (n[2]=0; n[2]<N[2]; ++(n[2]))
    for (n[1]=0; n[1]<N[1]; ++(n[1]))
      for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	for (int d=0; d<3; ++d) {
	  if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
	  else            nc[d] = n[d];
	  if (n[d]==Nh[d]) half[d] = 1;
	  else             half[d] = 0;
	}
	double scaled_rho = hypot(rhoR[nind], rhoI[nind])
	  /( (half[0]+1)*(half[1]+1)*(half[2]+1) );
	++nind;

	if (zero_vect(nc)) continue;
	if (scaled_rho < rho_min) continue;
        
	for (h[0]=0; h[0]<=half[0]; ++(h[0]))
	  for (h[1]=0; h[1]<=half[1]; ++(h[1]))
	    for (h[2]=0; h[2]<=half[2]; ++(h[2])) {
	      for (int d=0; d<3; ++d) nc_h[d] = nc[d] - h[d]*N[d];
	      Gabs = sqrt(magnsq(Gmetric, nc_h));
	      add_list(G_l, rho_l, Nl, log(Gabs), log(scaled_rho));
	    }
      }

  // now, fit to a power law...
  double G_sum=0, G2_sum=0, rho_sum=0, rhoG_sum=0;
  for (int i=0; i<Nl; ++i) {
    G_sum    += G_l[i];
    G2_sum   += G_l[i] * G_l[i];
    rho_sum  += rho_l[i];
    rhoG_sum += rho_l[i] * G_l[i];
  }
  filter_b = (Nl*rhoG_sum - G_sum*rho_sum)/(Nl*G2_sum - G_sum*G_sum);
  
  filter_a = 0;
  for (int i=0; i<Nl; ++i) {
    double guess_a =
      exp(rho_l[i] - filter_b*G_l[i]);
    if (guess_a > filter_a) filter_a = guess_a;
  }

  // garbage collection
  delete[] G_l;
  delete[] rho_l;
}



void apply_filter(double* rhoR, double* rhoI, 
		  int N[3], double latt[9],
		  double filter_a, double filter_b) 
{
  double Gmetric[9];
  {
    double temp[9];
    careful_inverse(latt, temp);	// temp = [a]^-1
    mult(temp, 2*M_PI, temp);		// temp = 2pi[a]^-1
    square(temp, Gmetric);		// Gmetric = (2pi)^2 [a]^-T[a]^-1
  }
    
  int n[3], nc[3];
  int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
  
  int nind=0;
  double Gabs, rho_magn, rho_max;
  
  for (n[2]=0; n[2]<N[2]; ++(n[2]))
    for (n[1]=0; n[1]<N[1]; ++(n[1]))
      for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	for (int d=0; d<3; ++d) {
	  if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
	  else            nc[d] = n[d];
	}
	rho_magn = hypot(rhoR[nind], rhoI[nind]);
	++nind;

	if (zero_vect(nc)) continue;
	Gabs = sqrt(magnsq(Gmetric, nc));
	rho_max = filter_a * exp(filter_b*log(Gabs));
	if (rho_magn > rho_max) {
	  rhoR[nind-1] *= rho_max/rho_magn;
	  rhoI[nind-1] *= rho_max/rho_magn;
	}
      }
}



// CROSS_LAMBDA[axis][0,1] : should have factor of lambda
// CROSS_LAMBDA[axis][2]   : remaining term, no factor of lambda
int CROSS_LAMBDA[3][3] = {
  {1,2, 5}, {5,1, 2}, {2,5, 1}};

// CROSS_LAMBDA[axis][0]   : should have factor of lambda^2
// CROSS_LAMBDA[axis][1,2] : remaining terms, no factor of lambda
int DIAG_LAMBDA[3][3] = {
  {0, 4,8}, {4, 8,0}, {8 ,0,4}};

// alpha = L1/(L1+L2) = 1/(1+L2/L1) = 1/(1 + 1/lambda)
int lattice_match(double latt0[9], double latt1[9], int axis, double &alpha) 
{
  double a[9], b[9];
  square(latt0, a);
  square(latt1, b);

  // No lambda terms:
  if ( (! dcomp(a[CROSS_LAMBDA[axis][2]],b[CROSS_LAMBDA[axis][2]])) ||
       (! dcomp(a[DIAG_LAMBDA[axis][1]],b[DIAG_LAMBDA[axis][1]])) ||
       (! dcomp(a[DIAG_LAMBDA[axis][2]],b[DIAG_LAMBDA[axis][2]])) )
    return 0;

  double lambda;
  lambda = sqrt(a[DIAG_LAMBDA[axis][0]]/b[DIAG_LAMBDA[axis][0]]);
  
  if ( (! dcomp(a[CROSS_LAMBDA[axis][0]],lambda*b[CROSS_LAMBDA[axis][0]])) ||
       (! dcomp(a[CROSS_LAMBDA[axis][1]],lambda*b[CROSS_LAMBDA[axis][1]])) )
    return 0;

  alpha = 1/(1+1/lambda);
  // if we made it this far, they agree
  return 1;
}
