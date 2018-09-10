/*
  Program: chgcar-disl.C
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: Using FFT, we'll take a CHGCAR file, read it's grid, and transform
           the file to include a screw dislocation.  Total kludge.

  Param.:  <CHGCAR> <u.1> <u.2> <angle>
           CHGCAR:  output from VASP; handles US, PAW, spin- and non-spin-polar
	   u.1,u.2: screw dislocation center in unit coord.
	   angle:   angle in DEGREES of dislocation cut relative to a1;
                    positive (or 0) means b = a3, negative means b = -a3

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   We read in everything, and then:
           1. Determine the displacement field for atoms, and output
	   2. Read the existing density grid, and compute a low-pass filter
	   3. Displace grid using cubic splines
	   4. FFT, apply a low-pass filter, and iFFT back--output
	   5. Output augmentation charges (if they exist)
	   6. repeat if spin-polarized

  Output:  New CHGCAR file for dislocation.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "3dfft.H"
#include "spline.H"
#include "chgcar.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

void extract_filter(double* rhoR, double* rhoI, 
		    int N[3], double latt[9],
		    double& filter_a, double& filter_b);

void apply_filter(double* rhoR, double* rhoI, 
		  int N[3], double latt[9],
		  double filter_a, double filter_b);

void apply_displacement(double* rho, int N[3],
			double D0, double D1, 
			double A0, double A1,
			double thetaD, int bsign);


inline double displace(const double& ux, const double& uy, const double& u0) 
{
  double theta =
    atan2(uy, ux) * (1./(2*M_PI)) - u0;
  if (theta < 0) return (theta+1);
  if (theta >= 1.) return (theta-1);
  return(theta);
}

inline double incell(const double& x) 
{
  double u = x - (int)(x);
  if (u<-0.5) return(u+1);
  if (u>=0.5) return(u-1);
  return(u);
}

inline double inunit(const double& x) 
{
  double u = x - (int)(x);
  if (u< 0) return(u+1);
  if (u>=1) return(u-1);
  return(u);
}

/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 4;
const char* ARGLIST = "<CHGCAR> <u.1> <u.2> <angle>";

const int NFLAGS = 0;
const char USERFLAGLIST[NFLAGS] = {}; // flag characters.

const char* ARGEXPL = 
"  CHGCAR:  output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  u.1,u.2: screw dislocation center in unit coord.\n\
  angle:   angle in DEGREES of dislocation cut relative to a1;\n\
           positive (or 0) means b = a3, negative means b = -a3";

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

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile;

  char* chgcar_name = args[0];
  double D0, D1, thetaD;
  int bsign;
  sscanf(args[1], "%lf", &D0);
  sscanf(args[2], "%lf", &D1);
  sscanf(args[3], "%lf", &thetaD);

  if (thetaD>=0) bsign=1;
  else           bsign=-1;
  
  thetaD *= 1./360.; // convert to radians, divided by 2pi
  for ( ; thetaD < 0; thetaD += 1.) ;
  for ( ; thetaD >= 1; thetaD -= 1.) ;
  
  
  infile = myopenr(chgcar_name);
  if (infile == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_name);
    exit(ERROR_NOFILE);
  }


  // ***************************** ANALYSIS **************************

  //++ comment line
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ a0
  fgets(dump, sizeof(dump), infile);
  printf("%s", dump);

  //++ lattice vectors
  double latt[9];
  for (i=0; i<3; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", latt+i, latt+3+i, latt+6+i);
    printf("%s", dump);
  }

  double amagn[9];
  square(latt, amagn);
  if ( (! zero(amagn[1])) || (! zero(amagn[2])) || (! zero(amagn[5])) ) {
    fprintf(stderr, "You're trying this with a non-rectangular lattice... but we're not ready for that.\n");
    exit(2);
  }

  double A0 = sqrt(amagn[0]), A1 = sqrt(amagn[4]);
  double u0, u1;

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
  double in_vect[3];
  for (i=0; i<Natoms; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", in_vect, in_vect+1, in_vect+2);
    u0 = incell(in_vect[0] - D0) * A0;
    u1 = incell(in_vect[1] - D1) * A1;
    in_vect[2] = inunit(in_vect[2] + bsign*displace(u0, u1, thetaD));
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
    
    // Now, all of our analysis.
    // 0. determine our low-pass filter values
    double* rhoR = new double[Ng];
    double* rhoI = new double[Ng];
    forward_3dfft(rho, N[0], N[1], N[2], rhoR, rhoI);
    double filter_a, filter_b;
    const double FILTER_MULT = 3.0;  // at some point, should make this a param...
    extract_filter(rhoR, rhoI, N, latt, filter_a, filter_b);
    if (VERBOSE) {
      fprintf(stderr, "# rho_filter(G)= %.8le |G|^ %.8lf\n", filter_a, filter_b);
      fprintf(stderr, "# mult_fact= %.8lf\n", FILTER_MULT);
    }
    filter_a *= FILTER_MULT;

    // 1. apply our displacement field
    apply_displacement(rho, N, D0, D1, A0, A1, thetaD, bsign);

    // 2. low-pass truncation filter
    forward_3dfft(rho, N[0], N[1], N[2], rhoR, rhoI);
    apply_filter(rhoR, rhoI, N, latt, filter_a, filter_b);
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



void apply_displacement(double* rho, int N[3],
			double D0, double D1, 
			double A0, double A1,
			double thetaD, int bsign) 
{
  int n0, n1, n2;
  double dN[3] = {1./N[0], 1./N[1], 1./N[2]};
  double u0, u1, du;

  int Np = N[2];
  double data[Np];
  cubic_type* spl = NULL;
  double dx = dN[2];

  int N0 = N[0], N1 = N[1], N2 = N[2], N01 = N[0]*N[1];
  for (n1=0; n1<N1; ++n1) 
    for (n0=0; n0<N0; ++n0) {
      u0 = incell((n0+0.5)*dN[0] - D0) * A0;
      u1 = incell((n1+0.5)*dN[1] - D1) * A1;
      du = bsign*displace(u0, u1, thetaD);
      int nind = n0+N0*n1;
      for (n2=0; n2<N2; ++n2)
	data[n2] = rho[nind + N01*n2];
      construct_spline(Np, dx, data, spl);
      for (n2=0; n2<N2; ++n2)
	rho[nind+N01*n2] = eval_spline(Np, dx, spl, n2*dN[2] - du);
    }
  delete[] spl;
}
