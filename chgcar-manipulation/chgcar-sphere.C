/*
  Program: chgcar-sphere.C
  Author:  D. Trinkle
  Date:    2007 July 25
  Purpose: We read in a bulk base CHGCAR (single atom unit cell), and attempt
           to construct a "spherical averaged" and smoothed version of the
	   charge density via linear interpolation

  Param.:  CHGCAR.bulk Ns [rcut] [smooth]
	   CHGCAR.bulk: output from VASP; handles US, PAW, spin- and non-spin-polar
           Ns:          number of "spline" points
	   rcut:        cutoff distance for spline; default = 3*(WS radius)
	   smooth:      weight factor for roughness penalty

  Flags:   MEMORY:  not used
	   VERBOSE: not used
	   TESTING: usual screen diahrea

  Algo.:   

  Output:  Spherical density evaluated at grid points
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
#include "chgcar.H"
#include <gsl/gsl_linalg.h>
  // #include "3dfft.H"
  // #include "linear.H"

//****************************** STRUCTURES ****************************


//***************************** SUBROUTINES ****************************
inline double cube_root (const double &x) { return exp(log(x)/3.); }

double nn(double latt[9]);

inline double Delta(const double &u) 
{
  if (u>1) return 0;
  if (u<-1) return 0;
  if (u>=0) return 1-u;
  return 1+u;
}


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "[-vt] CHGCAR.bulk Ns [rcut] [smooth]";

const char* ARGEXPL = 
"  CHGCAR.bulk: output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  Ns:          number of \"spline\" points\n\
  rcut:        cutoff distance for spline; default = 3*(WS radius)\n\
  smooth:      weight factor for roughness penalty\n\
\n\
  -v  VERBOSE\n\
  -t  TESTING\n";

const char* FILEEXPL =
"Both CHGCAR.bulk in VASP format";

int VERBOSE = 0;	// The infamous verbose flag.
int TESTING = 0;	// Extreme verbosity (testing purposes)
int ERROR = 0;		// Analysis: Error flag (for analysis purposes)

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);

  char ch;
  while ((ch = getopt(argc, argv, "vth")) != -1) {
    switch (ch) {
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
  FILE* infile_bulk;

  char* chgcar_bulk_name = argv[0];
  int Ns = strtol(argv[1], (char**)NULL, 0);
  double rcut = 0;
  if (argc > NUMARGS) rcut = strtod(argv[2], (char**)NULL);
  double smooth=0;
  if (argc > NUMARGS+1) smooth = strtod(argv[3], (char**)NULL);
  if (rcut < 0) {
    fprintf(stderr, "rcut=%g < 0?\n", rcut);
    exit(ERROR_BADFLAG);
  }
  
  infile_bulk = myopenr(chgcar_bulk_name);
  if (infile_bulk == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_bulk_name);
    exit(ERROR_NOFILE);
  }

  // ***************************** ANALYSIS **************************
  double a0;
  double latt_bulk[9], vol;
  int atomlist[100];
  int Natom=0;
  double **uatom = NULL;

  read_CHGCAR_header(infile_bulk, dump, a0, latt_bulk, atomlist, uatom);
  for (int i=0; atomlist[i]>=0; ++i) Natom += atomlist[i];
  for (int d=0; d<9; ++d) latt_bulk[d] *= a0;
  vol = det(latt_bulk);
  if ( (Natom != 1) ||
       (uatom[0][0] != 0) || (uatom[0][1] != 0) || (uatom[0][2] != 0) ) {
    fprintf(stderr, "Currently, we cannot handle starting from non-Bravais lattice bulk\n");
    exit(ERROR_BADFILE);
  }

  // cutoff:
  if (rcut == 0) 
    // Wigner-Seitz radius:
    rcut = 3.*cube_root( 3./(4.*M_PI) * vol );

  // Spline setup
  double dr=rcut/Ns, dr_1 = 1/dr;
  double *DD = new double[Ns*Ns], *RR = new double[Ns*Ns], *rhoD = new double[Ns];
  double *intdens = new double[Ns];
  double latt2[9];
  double nn_bulk = nn(latt_bulk);
  int Nmax = ceil(rcut/nn_bulk) + 1;

  for (int i=0; i<Ns; ++i) {
    double ri=(i*dr);
    intdens[i] = 2*M_PI/3*dr*(6*ri*ri + dr*dr);
    if (i==0) intdens[i] *= 0.5;
  }

  if (TESTING) {
    fprintf(stderr, "# Ns= %d\n", Ns);
    fprintf(stderr, "# rcut= %g  dr= %g\n", rcut, dr);
    fprintf(stderr, "# nn= %g\n", nn_bulk);
    fprintf(stderr, "# smooth= %g  Nmax= %d\n", smooth, Nmax);
  }

  // Read in grids, and analyze
  for (int ISPIN=0; ISPIN<2; ++ISPIN) {
    //++ grid values
    int N[3];
    //++ grid values
    fgets(dump, sizeof(dump), infile_bulk);
    sscanf(dump, "%d %d %d", N+0, N+1, N+2);
    int Ng=N[0]*N[1]*N[2];
    double* rho = new double[Ng]; // our soon-to-be charge density
    // set up our magnitude--designed so that we can run over the indices
    // without first scaling by 1/N
    {
      double temp[9];
      for (int i=0; i<3; ++i)
	for (int d=0; d<3; ++d) temp[3*d+i] = latt_bulk[3*d+i]/N[i];
      square(temp, latt2);
    }

    // ++ read in our bulk grid
    read_grid(infile_bulk, rho, Ng);
    if (feof(infile_bulk)) {
      myclose(infile_bulk);
      return ERROR;
    }

    // ++roughness
    for (int d=0; d<Ns*Ns; ++d) RR[d]=0;
    for (int i=0; i<Ns; ++i) {
      double ri=(i*dr);
      RR[i*(Ns+1)] = (2*dr*dr*dr*dr + 20*dr*dr*ri*ri + 15*ri*ri*ri*ri)/(15*dr);
      if (i==0) RR[0] *= 0.5;
      if (i<(Ns-1))
	RR[i*(Ns+1)+1] = 
	  -(9*dr*dr*dr*dr + 40*dr*dr*dr*ri + 70*dr*dr*ri*ri 
	    + 60*dr*ri*ri*ri + 30*ri*ri*ri*ri)/(60*dr);
      if (i>0)
	RR[i*(Ns+1)-1] = 
	  -(9*dr*dr*dr*dr - 40*dr*dr*dr*ri + 70*dr*dr*ri*ri 
	    - 60*dr*ri*ri*ri + 30*ri*ri*ri*ri)/(60*dr);
    }
    // --roughness

    // construct our fit matrices... nasty.
    for (int d=0; d<Ns*Ns; ++d) DD[d]=0;
    for (int d=0; d<Ns; ++d) rhoD[d]=0;

    int n[3], nind=0;
    for (n[2]=0; n[2]<N[2]; ++(n[2]))
      for (n[1]=0; n[1]<N[1]; ++(n[1]))
	for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	  double rho_value = rho[nind] / vol; // make into a density
	  
	  int n1[3];
	  for (n1[0]=n[0]-Nmax*N[0]; n1[0]<Nmax*N[0]; n1[0]+=N[0])
	    for (n1[1]=n[1]-Nmax*N[1]; n1[1]<Nmax*N[1]; n1[1]+=N[1])
	      for (n1[2]=n[2]-Nmax*N[2]; n1[2]<Nmax*N[2]; n1[2]+=N[2]) {
		double u1=sqrt(magnsq(latt2, n1))*dr_1;
		if (u1 >= Ns) continue;
		int i1[2]={floor(u1),floor(u1)+1};
		if (i1[0]>=0) rhoD[i1[0]] += rho_value*Delta(u1-i1[0]);
		if (i1[1]<Ns) rhoD[i1[1]] += rho_value*Delta(u1-i1[1]);

		int n2[3];
		for (n2[0]=n[0]-Nmax*N[0]; n2[0]<Nmax*N[0]; n2[0]+=N[0])
		  for (n2[1]=n[1]-Nmax*N[1]; n2[1]<Nmax*N[1]; n2[1]+=N[1])
		    for (n2[2]=n[2]-Nmax*N[2]; n2[2]<Nmax*N[2]; n2[2]+=N[2]) {
		      double u2=sqrt(magnsq(latt2, n2))*dr_1;
		      if (u2 >= Ns) continue;
		      int i2[2]={floor(u2),floor(u2)+1};
		      for (int d1=0; d1<2; ++d1) 
			if ( (i1[d1]>=0) && (i1[d1]<Ns) ) 
			  for (int d2=0; d2<2; ++d2) 
			    if ( (i2[d2]>=0) && (i2[d2]<Ns) ) 
			      DD[i1[d1]*Ns+i2[d2]] += Delta(u1-i1[d1])*Delta(u2-i2[d2]);
		    }
	      }
	  ++nind;
	}
    if (TESTING) {
      fprintf(stderr, "# DD (onsite):\n");
      for (int i=0; i<Ns; ++i) {
	fprintf(stderr, "# DD[%g]= %g", i*dr, DD[i*(Ns+1)]);
	fprintf(stderr, " DD+[%g]= %g", i*dr, DD[i*(Ns+1)+1]);
	fprintf(stderr, " DD-[%g]= %g", i*dr, DD[i*(Ns+1)-1]);
	fprintf(stderr, "\n");
      }
      fprintf(stderr, "# RR (onsite):\n");
      for (int i=0; i<Ns; ++i) {
	fprintf(stderr, "# RR[%g]= %g", i*dr, RR[i*(Ns+1)]);
	fprintf(stderr, " RR+[%g]= %g", i*dr, RR[i*(Ns+1)+1]);
	fprintf(stderr, " RR-[%g]= %g", i*dr, RR[i*(Ns+1)-1]);
	fprintf(stderr, "\n");
      }
      
    }
    // now, solve...
    {
      // 1. solve with no roughness penalty first:
      double *FF = new double[Ns*Ns];
      for (int d=0; d<Ns*Ns; ++d) FF[d] = DD[d];
      gsl_matrix_view m = gsl_matrix_view_array (FF, Ns, Ns);
      gsl_vector_view b = gsl_vector_view_array (rhoD, Ns);
      gsl_vector *x = gsl_vector_alloc (Ns);
      int s;
      gsl_permutation * p = gsl_permutation_alloc (Ns);
      gsl_linalg_LU_decomp (&m.matrix, p, &s);
      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

      // 2. calculate chi^2 and initial roughness:
      double chi2=0, r2=0;
      for (int ng=0; ng<Ng; ++ng) chi2 += 0.5*rho[ng]*rho[ng]/(vol*vol);
      for (int i=0; i<Ns; ++i)
	chi2 += -0.5*rhoD[i]*gsl_vector_get(x,i);

      for (int i=0; i<Ns; ++i)
	for (int j=0; j<Ns; ++j)
	  r2 += gsl_vector_get(x,i)*RR[i*Ns+j]*gsl_vector_get(x,j);
      if (TESTING) {
	fprintf(stderr, "# chi2= %g  roughness= %g\n", chi2, r2);
      }
      // 3. resolve with roughness penalty:
      double alpha = smooth*chi2/r2;
      for (int d=0; d<Ns*Ns; ++d) FF[d] = DD[d] + alpha*RR[d];
      gsl_linalg_LU_decomp (&m.matrix, p, &s);
      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
      // 4. scale to get total density correct
      double charge=0, chargefit=0;
      for (int ng=0; ng<Ng; ++ng) charge += rho[ng];
      charge /= Ng;
      for (int i=0; i<Ns; ++i) chargefit += intdens[i] * gsl_vector_get(x, i);
      chi2=0, r2=0;
      for (int ng=0; ng<Ng; ++ng) chi2 += 0.5*rho[ng]*rho[ng]/(vol*vol);
      for (int i=0; i<Ns; ++i)
	chi2 += -0.5*rhoD[i]*gsl_vector_get(x,i);

      for (int i=0; i<Ns; ++i)
	for (int j=0; j<Ns; ++j)
	  r2 += gsl_vector_get(x,i)*RR[i*Ns+j]*gsl_vector_get(x,j);
      if (TESTING) {
	fprintf(stderr, "# chi2= %g  roughness= %g\n", chi2, r2);
      }
      if (TESTING) {
	fprintf(stderr, "# charge= %g  chargefit= %g\n", charge, chargefit);
      }
      // 5. output
      printf("%d # number of grid points\n", Ns);
      for (int ns=0; ns<Ns; ++ns)
	printf("%g %g\n", ns*dr, gsl_vector_get(x, ns) * charge/chargefit);
      printf("%g %g\n", rcut, 0.0);
      gsl_permutation_free (p);
      gsl_vector_free (x);
    }
    
    //++ (augmentation charges?)
    fgets(dump, sizeof(dump), infile_bulk);
    if (dump[0] == 'a') {
      // augmentation charges...
      int aug_grid;
      sscanf(dump, "%*s %*s %*d %d", &aug_grid);
      double* aug = new double[aug_grid];
      read_grid(infile_bulk, aug, aug_grid);
        
      printf("augmentation occupancies %3d %3d\n", 1, aug_grid);
      write_grid(stdout, aug, aug_grid);
      fgets(dump, sizeof(dump), infile_bulk);
    }
    
    if (feof(infile_bulk)) {
      myclose(infile_bulk);
      return ERROR;
    }
  
    //++ we're spin-polarized, so we have to read the magmom grid, output, and go back.
    double magmom;
    // ... parse the line we've already read.
    magmom=strtod(dump, (char**)NULL);
    // printf magmom...
  }
  
  myclose(infile_bulk);

  delete[] DD;
  delete[] RR;
  delete[] rhoD;
  delete[] intdens;

  return ERROR;
}


// calculates (and returns) the nearest neighbor distance
double nn(double latt[9]) 
{
  double latt2[9];
  square(latt, latt2);
  double rmin = latt2[0];
  int n[3];
  for (n[0]=-2; n[0]<=2; ++n[0])
    for (n[1]=-2; n[1]<=2; ++n[1])
      for (n[2]=-2; n[2]<=2; ++n[2]) {
	double r2=magnsq(latt2, n);
	if ( (r2>0) && (r2 < rmin) ) rmin=r2;
      }
  return sqrt(rmin);
}
