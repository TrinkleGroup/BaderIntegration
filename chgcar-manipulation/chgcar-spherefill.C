/*
  Program: chgcar-spherefill.C
  Author:  D. Trinkle
  Date:    2007 July 27
  Purpose: We read in a spherically averaged charge density, and the
           beginning of a new CHGCAR file containing atoms in a new arrangement,
	   with a specified grid.  We construct a new CHGCAR file for the new
	   arrangement, by superimposing our spherically averaged charge densities.
	   
	   This charge density is intended to be an initial guess for DFT
	   calculations of a defect.

  Param.:  CHGCAR.sphere CHGCAR.new
	   CHGCAR.sphere: output from chgcar-sphere (fit from bulk CHGCAR)
           CHGCAR.new:    header portion of CHGCAR, including grid

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   

  Output:  New CHGCAR file, now with density.
*/

//************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include "io.H"    // All of our "read in file", etc.
#include "dcomp.H"
#include "matrix.H"
#include "3dfft.H"
#include "linear.H"
#include "chgcar.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************

inline double Delta(const double &u) 
{
  if (u>1) return 0;
  if (u<-1) return 0;
  if (u>=0) return 1-u;
  return 1+u;
}

// r:    point to evaluate
// Ns:   number of spline points
// dr_1: 1/dr (inverse spacing)
// f_i:  array of spline points
inline double f_r (const double &r, 
		   const int &Ns, const double &dr_1,
		   double f_i[]) 
{
  double u = r*dr_1;
  int i = (int)u;
  if (i<0) return f_i[0];
  if (i>=Ns) return 0;
  double f = f_i[i]*Delta(u-i);
  if (i<(Ns-1)) f += f_i[i+1]*Delta(u-i-1);
  return f;
}

void populate_grid(double* rho, int N[3], double latt[9], 
		   int Natom, double** uatom,
		   const int &Ns, const double &dr, double rho_i[]);

inline double cube_root (double x) { return exp(log(x)/3.); }


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "[-krvt] CHGCAR.sphere CHGCAR.new";

const char* ARGEXPL = 
"  CHGCAR.sphere: output from chgcar-sphere (fit from bulk CHGCAR)\n\
  CHGCAR.new:    header portion of CHGCAR, including grid\n\
\n\
  -k  keep charge density; no \"resetting\" to correct amount\n\
  -v  VERBOSE\n\
  -t  TESTING\n";

const char* FILEEXPL =
"CHGCAR.new is in VASP format; CHGCAR.new only requires\n\
the atom positions and the new grid dimensions.\n\
==== CHCGCAR.sphere ====\n\
[number of spline points]\n\
r_0 rho_0\n\
...\n\
r_(N-1) rho_(N-1)\n\
==== CHGCAR.sphere ====\n";

int VERBOSE = 0;	// The infamous verbose flag.
int TESTING = 0;	// Extreme verbosity (testing purposes)
int ERROR = 0;		// Analysis: Error flag (for analysis purposes)

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);
  int RESCALE = 1;	// rescale total charge

  char ch;
  while ((ch = getopt(argc, argv, "kvth")) != -1) {
    switch (ch) {
    case 'k':
      RESCALE = 0;
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

  if (argc < NUMARGS) ERROR = 2;
  // All hell broken loose yet?
  if (ERROR != 0) {
    fprintf(stderr, "%s %s\n", progname, ARGLIST);
    fprintf(stderr, "%s\n", ARGEXPL);
    fprintf(stderr, "Input file format:\n%s\n", FILEEXPL);
    exit(ERROR);
  }

  // ****************************** INPUT ****************************
  char dump[512];
  FILE* infile_sphere;
  FILE* infile_new;

  char* chgcar_sphere_name = argv[0];
  char* chgcar_new_name = argv[1];

  infile_sphere = myopenr(chgcar_sphere_name);
  if (infile_sphere == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_sphere_name);
    exit(ERROR_NOFILE);
  }

  infile_new = myopenr(chgcar_new_name);
  if (infile_new == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_new_name);
    exit(ERROR_NOFILE);
  }

  if (infile_sphere == infile_new) {
    fprintf(stderr, "Cannot have both CHGCAR.sphere and CHGCAR.new be stdin.\n");
    exit(ERROR_NOFILE);
  }
  
  // ***************************** ANALYSIS **************************
  double a0;
  double latt_new[9];
  int atomlist[100];
  int Natom_new=0;
  double **uatom_new = NULL;

  read_CHGCAR_header(infile_new, dump, a0, latt_new, atomlist, uatom_new);
  write_CHGCAR_header(stdout, dump, a0, latt_new, atomlist, uatom_new);
  for (int i=0; atomlist[i]>=0; ++i) Natom_new += atomlist[i];
  for (int d=0; d<9; ++d) latt_new[d] *= a0;

  //++ grid values
  int N[3];
  nextnoncomment(infile_new, dump, sizeof(dump));
  sscanf(dump, "%d %d %d", N+0, N+1, N+2);
  int Ng=N[0]*N[1]*N[2];
  double* rho = new double[Ng]; // our soon-to-be charge density
  myclose(infile_new);
  
  for (int ISPIN=0; ISPIN<2; ++ISPIN) {
    //++ grid values
    nextnoncomment(infile_sphere, dump, sizeof(dump));
    int Ns = strtol(dump, (char**)NULL, 0);
    double *r_i=new double[Ns+1], *rho_i=new double[Ns+1];
    for (int n=0; n<=Ns; ++n) {
      nextnoncomment(infile_sphere, dump, sizeof(dump));
      sscanf(dump, "%lf %lf", r_i+n, rho_i+n);
    }
    double dr = (r_i[Ns]-r_i[0])/Ns;
    double charge=0;
    for (int n=0; n<Ns; ++n) {
      if (n==0) charge += M_PI*dr*dr*dr/3 * rho_i[0];
      else 
	charge += 2*M_PI/3 * dr*(6*r_i[n]*r_i[n] + dr*dr) * rho_i[n];
    }

    // Now, all of our analysis.
    // 0. populate the 
    populate_grid(rho, N, latt_new, Natom_new, uatom_new, Ns, dr, rho_i);

    // garbage collection
    delete[] r_i;
    delete[] rho_i;

    // quick sanity check against "negative" charge:
    if (ISPIN == 0) {
      double basescale = (double)Natom_new * (double)Ng;
      for (int ng=0; ng<Ng; ++ng)
	if (rho[ng]<0) rho[ng]=0;
      if (RESCALE) {
	double cell_new=0;
	for (int ng=0; ng<Ng; ++ng) cell_new += rho[ng];
	if (TESTING) {
	  fprintf(stderr, "# sum rho=  %g x %d x %d\n", 
		  cell_new/basescale, Natom_new, Ng);
	  fprintf(stderr, "# should be= %g\n", charge);
	}
	double scale = charge * basescale/ cell_new;
	for (int ng=0; ng<Ng; ++ng) rho[ng] *= scale;
      }
    }

    // *** OUTPUT ***
    printf(" %4d %4d %4d\n", N[0], N[1], N[2]);
    write_grid(stdout, rho, Ng);
    
    if (feof(infile_sphere)) {
      myclose(infile_sphere);
      return ERROR;
    }
    
    //++ (augmentation charges?)
    nextnoncomment(infile_sphere, dump, sizeof(dump));
    if (dump[0] == 'a') {
      // augmentation charges...
      int aug_grid;
      sscanf(dump, "%*s %*s %*d %d", &aug_grid);
      double* aug = new double[aug_grid];
      read_grid(infile_sphere, aug, aug_grid);
        
      for (int n=0; n<Natom_new; ++n) {
	printf("augmentation occupancies %3d %3d\n", n+1, aug_grid);
	write_grid(stdout, aug, aug_grid);
      }
      nextnoncomment(infile_sphere, dump, sizeof(dump));
    }
    
    if (feof(infile_sphere)) {
      myclose(infile_sphere);
      return ERROR;
    }
  
    //++ we're spin-polarized, so we have to read the magmom grid, output, and go back.
    double magmom, *magmom_mult = new double[Natom_new];
    // ... parse the line we've already read.
    magmom=strtod(dump, (char**)NULL);
    for (int n=0; n<Natom_new; ++n)
        magmom_mult[n] = magmom;
    write_grid(stdout, magmom_mult, Natom_new);
    delete[] magmom_mult;

  }
  
  myclose(infile_sphere);

  return ERROR;
}

// Brute force approach to determine the possible shifts that need to be considered
void calc_shifts(double latt[9], int Nshift[3], const double &rcut) 
{
  double metric[9];
  square(latt, metric);
  for (int d=0; d<3; ++d) Nshift[d]=0;
  double rc2=rcut*rcut;
  int shift[3];
  int Nmax[3]={2*ceil(sqrt(metric[0]/rc2)),
	       2*ceil(sqrt(metric[4]/rc2)),
	       2*ceil(sqrt(metric[8]/rc2))};
  for (shift[0]=-Nmax[0]; shift[0]<=Nmax[0]; ++shift[0])
    for (shift[1]=-Nmax[1]; shift[1]<=Nmax[1]; ++shift[1])
      for (shift[2]=-Nmax[2]; shift[2]<=Nmax[2]; ++shift[2]) {
	double rmagn2 = magnsq(metric, shift);
	if (rmagn2 < rc2) 
	  for (int d=0; d<3; ++d)
	    if (Nshift[d] < abs(shift[d])) Nshift[d] = abs(shift[d]);
      }
  for (int d=0; d<3; ++d) ++(Nshift[d]);
}


void populate_grid(double* rho, int N[3], double latt[9], 
		   int Natom, double** uatom,
		   const int &Ns, const double &dr, double rho_i[])
{
  double metric[9];
  square(latt, metric);

  int N0=N[0], N1=N[1], N2=N[2];
  double iN0=1./(double)N0, iN1=1./(double)N1, iN2=1./(double)N2;

  double rcut2 = (Ns*Ns*dr*dr);
  double dr_1 = 1./dr;
  
  int Nshift[3];
  calc_shifts(latt, Nshift, sqrt(rcut2));

  double u[3];
  double chargescale = fabs(det(latt));
  int ng=0, NG=N0*N1*N2, NGprint=NG/64;
  if (VERBOSE) fprintf(stderr, "# populating...\n");
  for (int n2=0; n2<N2; ++n2) {
    u[2]=n2*iN2;
    for (int n1=0; n1<N1; ++n1) {
      u[1]=n1*iN1;
      for (int n0=0; n0<N0; ++n0) {
	u[0]=n0*iN0;
	int shift[3];
	for (int n=0; n<Natom; ++n) {
	  double ud0[3] = {u[0]-uatom[n][0],u[1]-uatom[n][1],u[2]-uatom[n][2]};
	  for (shift[0]=-Nshift[0]; shift[0]<=Nshift[0]; ++shift[0])
	    for (shift[1]=-Nshift[1]; shift[1]<=Nshift[1]; ++shift[1])
	      for (shift[2]=-Nshift[2]; shift[2]<=Nshift[2]; ++shift[2]) {
		double ud[3] = {ud0[0]+shift[0],ud0[1]+shift[1],ud0[2]+shift[2]};
		double rmagn2 = magnsq(metric, ud);
		if (rmagn2 < rcut2) 
		  // add in charge density here
		  rho[ng] += chargescale * f_r(sqrt(rmagn2), Ns, dr_1, rho_i);
	      }
	}
	++ng;
	if (VERBOSE)
	  if (ng%NGprint == 0) fprintf(stderr, "#");
      }
    }
  }
  if (VERBOSE) fprintf(stderr, "\n");
}
