/*
  Program: chgcar-voronoi.C
  Author:  D. Trinkle
  Date:    2007 July 5
  Purpose: We read in a bulk base CHGCAR (single atom unit cell), and the
           beginning of a new CHGCAR file containing atoms in a new arrangement,
	   with a specified grid.  We construct a new CHGCAR file for the new
	   arrangement, by identifying each point on the new CHGCAR grid
	   with one atom a la the Voronoi tesselation, and then determining
	   that point's charge density from the bulk base CHGCAR file.  If
	   the grid point is sufficiently far from any atom (i.e., in "vacuum") we
	   assume the point has zero charge density.
	   
	   The final charge density is smoothed via a low-pass filter after a
	   FT, where the number of electrons is reset to the correct value, and
	   then IFTed back.

	   This charge density is intended to be an initial guess for DFT
	   calculations of a defect.

  Param.:  CHGCAR.bulk CHGCAR.new [rcut]
	   CHGCAR.bulk: output from VASP; handles US, PAW, spin- and non-spin-polar
           CHGCAR.new:  header portion of CHGCAR, including grid
	   rcut:        cutoff distance for being "vacuum"; default = WS radius

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

void populate_grid(double* rho, int N[3], double latt[9], 
		   int Natom, double** uatom,
		   double* rhob, int Nb[3], double latt_bulk[9],
		   double rcut, double rsmall);

double min_nn(double latt[9], int N, double** u);


void extract_filter(double* rhoR, double* rhoI, 
		    int N[3], double latt[9],
		    double& filter_a, double& filter_b);

void apply_filter(double* rhoR, double* rhoI, 
		  int N[3], double latt[9],
		  double filter_a, double filter_b);


void extract_filter_gaussian(double* rhoR, double* rhoI, 
			     int N[3], double latt[9],
			     double& l0);

void apply_filter_gaussian(double* rhoR, double* rhoI, 
			   int N[3], double latt[9],
			   double l0);

inline double inunit(const double& x) 
{
  double u = x - (int)(x);
  if (u< 0) return(u+1);
  if (u>=1) return(u-1);
  return(u);
}

inline double cube_root (double x) { return exp(log(x)/3.); }


/*================================= main ==================================*/

// Arguments first, then flags, then explanation.
const int NUMARGS = 2;
const char* ARGLIST = "[-rvt] CHGCAR.bulk CHGCAR.new [rcut]";

const char* ARGEXPL = 
"  CHGCAR.bulk: output from VASP; handles US, PAW, spin- and non-spin-polar\n\
  CHGCAR.new:  header portion of CHGCAR, including grid\n\
  rcut:        cutoff distance for being \"vacuum\"; default = WS radius\n\
\n\
  -k  keep charge density; no \"resetting\" to correct amount\n\
  -r  \"raw\" density -- no smoothing\n\
  -v  VERBOSE\n\
  -t  TESTING\n\
\n\
** WARNING: CHGCAR.bulk and CHGCAR.new should _not_ both be stdin";

const char* FILEEXPL =
"Both CHGCAR.bulk and CHGCAR.new are in VASP format; CHGCAR.new only requires\n\
the atom positions and the new grid dimensions, while CHGCAR.bulk must be a\n\
complete CHGCAR file, output from VASP to be meaningful.";

int VERBOSE = 0;	// The infamous verbose flag.
int TESTING = 0;	// Extreme verbosity (testing purposes)
int ERROR = 0;		// Analysis: Error flag (for analysis purposes)

int main ( int argc, char **argv ) 
{
  // ************************** INITIALIZATION ***********************
  char* progname = basename(argv[0]);
  int RESCALE = 1;	// rescale total charge
  int SMOOTH = 1;	// smooth total charge

  char ch;
  while ((ch = getopt(argc, argv, "krvth")) != -1) {
    switch (ch) {
    case 'k':
      RESCALE = 0;
      break;
    case 'r':
      SMOOTH = 0;
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
  FILE* infile_bulk;
  FILE* infile_new;

  char* chgcar_bulk_name = argv[0];
  char* chgcar_new_name = argv[1];
  double rcut = 0;
  if (argc > NUMARGS) rcut = strtod(argv[2], (char**)NULL);
  if (rcut < 0) {
    fprintf(stderr, "rcut=%g < 0?\n", rcut);
    exit(ERROR_BADFLAG);
  }
  
  infile_bulk = myopenr(chgcar_bulk_name);
  if (infile_bulk == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_bulk_name);
    exit(ERROR_NOFILE);
  }

  infile_new = myopenr(chgcar_new_name);
  if (infile_new == NULL) {
    fprintf(stderr, "Couldn't open %s for reading.\n", chgcar_new_name);
    exit(ERROR_NOFILE);
  }

  if (infile_bulk == infile_new) {
    fprintf(stderr, "Cannot have both CHGCAR.bulk and CHGCAR.new be stdin.\n");
    exit(ERROR_NOFILE);
  }
  
  // ***************************** ANALYSIS **************************
  double a0;
  double latt_bulk[9], latt_new[9];
  int atomlist[100];
  int Natom_bulk=0, Natom_new=0;
  double **uatom_bulk = NULL, **uatom_new = NULL;
  double rsmall;

  read_CHGCAR_header(infile_bulk, dump, a0, latt_bulk, atomlist, uatom_bulk);
  for (int i=0; atomlist[i]>=0; ++i) Natom_bulk += atomlist[i];
  for (int d=0; d<9; ++d) latt_bulk[d] *= a0;
  if ( (Natom_bulk != 1) ||
       (uatom_bulk[0][0] != 0) || (uatom_bulk[0][1] != 0) || (uatom_bulk[0][2] != 0) ) {
    fprintf(stderr, "Currently, we cannot handle starting from non-Bravais lattice bulk\n");
    exit(ERROR_BADFILE);
  }
  read_CHGCAR_header(infile_new, dump, a0, latt_new, atomlist, uatom_new);
  write_CHGCAR_header(stdout, dump, a0, latt_new, atomlist, uatom_new);
  for (int i=0; atomlist[i]>=0; ++i) Natom_new += atomlist[i];
  for (int d=0; d<9; ++d) latt_new[d] *= a0;

  // cutoff:
  if (rcut == 0) 
    // Wigner-Seitz radius:
    rcut = cube_root( 3./(4.*M_PI) *det(latt_bulk) );
  // minimum distance:
  rsmall = 0.5*min_nn(latt_new, Natom_new, uatom_new);

  //++ grid values
  int N[3];
  fgets(dump, sizeof(dump), infile_new);
  sscanf(dump, "%d %d %d", N+0, N+1, N+2);
  int Ng=N[0]*N[1]*N[2];
  double* rho = new double[Ng]; // our soon-to-be charge density
  myclose(infile_new);
  
  for (int ISPIN=0; ISPIN<2; ++ISPIN) {
    //++ grid values
    int Nb[3];
    fgets(dump, sizeof(dump), infile_bulk);
    sscanf(dump, "%d %d %d", Nb+0, Nb+1, Nb+2);
    
    int Nbg=Nb[0]*Nb[1]*Nb[2];
    double* rhob = new double[Nbg];

    // ++ read in our bulk grid
    read_grid(infile_bulk, rhob, Nbg);

    // Now, all of our analysis.
    // 0. populate the 
    populate_grid(rho, N, latt_new, Natom_new, uatom_new, rhob, Nb, latt_bulk,
		  rcut, rsmall);

    // 1. determine our low-pass filter values (from bulk)
    if (TESTING) {
      double cell=0;
      for (int ng=0; ng<Nbg; ++ng) cell += rhob[ng];
      fprintf(stderr, "# sum rhob= %g x %d x %d\n", 
	      cell/(Natom_bulk*Nbg), Natom_bulk, Nbg);
      int negcount=0;
      for (int ng=0; ng<Nbg; ++ng) if (rhob[ng]<0) ++negcount;
      fprintf(stderr, "# negcount= %d/%d = %g\n", 
	      negcount, Nbg, (double)negcount/(double)Nbg);      
    }
    //    double filter_a, filter_b;
    double l0;
    double cell_charge;
    const double FILTER_MULT = 3.0;  // at some point, should make this a param...
    {
      double* rhobR = new double[Nbg];
      double* rhobI = new double[Nbg];
      forward_3dfft(rhob, Nb[0], Nb[1], Nb[2], rhobR, rhobI);
      cell_charge = rhobR[0]/Nbg; // need this for later
      // extract_filter_gaussian(rhobR, rhobI, Nb, latt_bulk, filter_a, filter_b);
      extract_filter_gaussian(rhobR, rhobI, Nb, latt_bulk, l0);
      if (VERBOSE) {
	fprintf(stderr, "# rsmall= %g\n", rsmall);
	fprintf(stderr, "# rcut=   %g\n", rcut);
	// fprintf(stderr, "# rho_filter(G)= %.8le |G|^ %.8lf\n", filter_a, filter_b);
	fprintf(stderr, "# rho l0= %g\n", l0);
	fprintf(stderr, "# mult_fact= %.8lf\n", FILTER_MULT);
	fprintf(stderr, "# scale_fact= %d\n", Natom_new);
      }
      // filter_a *= FILTER_MULT * Natom_new; // also scale by number of atoms
      // garbage collection
      delete[] rhob;
      delete[] rhobR;
      delete[] rhobI;
    }
    
    // 2. low-pass truncation filter
    double* rhoR = new double[Ng];
    double* rhoI = new double[Ng];
    double basescale = (double)Ng * (double)Natom_new; // needed due to size
    if (TESTING) {
      double cell_new=0;
      for (int ng=0; ng<Ng; ++ng) cell_new += rho[ng];
      fprintf(stderr, "# sum rho=  %g x %d x %d\n", 
	      cell_new/basescale, Natom_new, Ng);
      int negcount=0;
      for (int ng=0; ng<Ng; ++ng) if (rho[ng]<0) ++negcount;
      fprintf(stderr, "# negcount= %d/%d = %g\n", 
	      negcount, Ng, (double)negcount/(double)Ng);      
    }
    forward_3dfft(rho, N[0], N[1], N[2], rhoR, rhoI);
    if (VERBOSE) {
      fprintf(stderr, "# rhoR[G=0]= %g x %d x %d\n", 
	      rhoR[0]/basescale, Natom_new, Ng);
      fprintf(stderr, "# should be= %g\n", cell_charge);
      fprintf(stderr, "# RESCALE= %d\n", RESCALE);
      fprintf(stderr, "# SMOOTH= %d\n", SMOOTH);
    }
    // scale values to correct charge / spin:
    if (RESCALE) {
      double scale = cell_charge * basescale/ rhoR[0];
      for (int ng=0; ng<Ng; ++ng) {
	rhoR[ng] *= scale;
	rhoI[ng] *= scale;
      }
    }
    // smooth it out:
    // if (SMOOTH) apply_filter(rhoR, rhoI, N, latt_new, filter_a, filter_b);
    if (SMOOTH) apply_filter_gaussian(rhoR, rhoI, N, latt_new, l0/FILTER_MULT);
    inverse_3dfft(rhoR, rhoI, N[0], N[1], N[2], rho);
    // garbage collection
    delete[] rhoR;
    delete[] rhoI;

    // quick sanity check against "negative" charge:
    if (ISPIN == 0) {
      for (int ng=0; ng<Ng; ++ng)
	if (rho[ng]<0) rho[ng]=0;
      if (RESCALE) {
	double cell_new=0;
	for (int ng=0; ng<Ng; ++ng) cell_new += rho[ng];
	if (TESTING) {
	  fprintf(stderr, "# sum rho=  %g x %d x %d\n", 
		  cell_new/basescale, Natom_new, Ng);
	}
	double scale = cell_charge * basescale/ cell_new;
	for (int ng=0; ng<Ng; ++ng) rho[ng] *= scale;
      }
    }

    // *** OUTPUT ***
    printf(" %4d %4d %4d\n", N[0], N[1], N[2]);
    write_grid(stdout, rho, Ng);
    
    if (feof(infile_bulk)) {
      myclose(infile_bulk);
      return ERROR;
    }
    
    //++ (augmentation charges?)
    fgets(dump, sizeof(dump), infile_bulk);
    if (dump[0] == 'a') {
      // augmentation charges...
      int aug_grid;
      sscanf(dump, "%*s %*s %*d %d", &aug_grid);
      double* aug = new double[aug_grid];
      read_grid(infile_bulk, aug, aug_grid);
        
      for (int n=0; n<Natom_new; ++n) {
	printf("augmentation occupancies %3d %3d\n", n+1, aug_grid);
	write_grid(stdout, aug, aug_grid);
      }
      fgets(dump, sizeof(dump), infile_bulk);
    }
    
    if (feof(infile_bulk)) {
      myclose(infile_bulk);
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
  
  myclose(infile_bulk);

  return ERROR;
}

// 3^3 * N*(N+1)/2 algorithm for finding shortest NN distance
double min_nn(double latt[9], int N, double** u) 
{
  double metric[9];
  square(latt, metric);
  double rmin2=metric[0];
  for (int m=0; m<N; ++m)
    for (int n=m; n<N; ++n) {
      double ud0[3] = {u[n][0]-u[m][0],u[n][1]-u[m][1],u[n][2]-u[m][2]}, ud[3];
      int shift[3];
      for (shift[0]=-1; shift[0]<=1; ++shift[0])
	for (shift[1]=-1; shift[1]<=1; ++shift[1])
	  for (shift[2]=-1; shift[2]<=1; ++shift[2]) {
	    for (int d=0; d<3; ++d) ud[d] = ud0[d] + shift[d];
	    double rmagn2 = magnsq(metric, ud);
	    if (rmagn2 == 0) continue;
	    if (rmagn2 < rmin2) rmin2 = rmagn2;
	  }
    }
  return sqrt(rmin2);
}


void populate_grid(double* rho, int N[3], double latt[9], 
		   int Natom, double** uatom,
		   double* rhob, int Nb[3], double latt_bulk[9],
		   double rcut, double rsmall)
{
  double metric[9];
  square(latt, metric);

  double invlatt[9]; // = [a]^-1 [A] ; turns coord. in supercell into bulk position
  { double invlatt_bulk[9];
    careful_inverse(latt_bulk, invlatt_bulk);
    mult(invlatt_bulk, latt, invlatt);
  }
  
  int N0=N[0], N1=N[1], N2=N[2];
  double iN0=1./(double)N0, iN1=1./(double)N1, iN2=1./(double)N2;

  double rsmall2=rsmall*rsmall;
  
  double u[3];
  double chargescale = fabs(det(latt)/det(latt_bulk)); // note: _not_ Natom
  int ng=0, NG=N0*N1*N2, NGprint=NG/64;
  if (VERBOSE) fprintf(stderr, "# populating...\n");
  for (int n2=0; n2<N2; ++n2) {
    u[2]=n2*iN2;
    for (int n1=0; n1<N1; ++n1) {
      u[1]=n1*iN1;
      for (int n0=0; n0<N0; ++n0) {
	u[0]=n0*iN0;
	double rmin2=rcut*rcut, umin[3];
	int atom=-1;
	int shift[3];
	for (int n=0; (n<Natom) && (rmin2 > rsmall2); ++n) {
	  double ud0[3] = {u[0]-uatom[n][0],u[1]-uatom[n][1],u[2]-uatom[n][2]},
	    ud[3];
	  for (shift[0]=-1; shift[0]<=1; ++shift[0])
	    for (shift[1]=-1; shift[1]<=1; ++shift[1])
	      for (shift[2]=-1; shift[2]<=1; ++shift[2]) {
		for (int d=0; d<3; ++d) ud[d] = ud0[d] + shift[d];
		double rmagn2 = magnsq(metric, ud);
		if (rmagn2 < rmin2) {
		  rmin2 = rmagn2;
		  atom = n;
		  for (int d=0; d<3; ++d) umin[d] = ud[d];
		}
	      }
	}
	if (atom < 0) 
	  rho[ng] = 0;
	else {
	  // linear interpolation to find -- have to mult by Natom, too.
	  double ubulk[3];
	  mult_vect(invlatt, umin, ubulk);
	  for (int d=0; d<3; ++d) ubulk[d] = inunit(ubulk[d])*Nb[d];
	  rho[ng] = chargescale * linear(rhob, Nb, ubulk);
	}
	++ng;
	if (VERBOSE)
	  if (ng%NGprint == 0) fprintf(stderr, "#");
      }
    }
  }
  if (VERBOSE) fprintf(stderr, "\n");
}


// *************************** POWERLAW FILTERS ************************

void add_list (double* G_l, double* rho_l, int& Nl, double G, double rho) 
{
  // I don't know why this was here... I've removed it
  //  if (zero(rho)) return;

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

const double RHO_MIN_SCALE = 1e-4;

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
  
  int n[3], nc[3];
  int Nh[3] = {N[0]/2, N[1]/2, N[2]/2};
  
  int nind=0;
  int h[3], half[3];
  double Gabs;
  int MEMORY = (Nh[0]+1)*(Nh[1]+1)*(Nh[2]+1); // hardcoded, I know...
  double* G_l = new double[MEMORY];
  double* rho_l = new double[MEMORY];
  int Nl = 0;

  double rho_min = rhoR[0] * RHO_MIN_SCALE;
  
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
	      double nc_h[3];
	      for (int d=0; d<3; ++d) nc_h[d] = nc[d] - h[d]*N[d];
	      Gabs = sqrt(magnsq(Gmetric, nc_h));
	      add_list(G_l, rho_l, Nl, log(Gabs), log(scaled_rho));
	      if (Nl >= MEMORY) {
		fprintf(stderr, "** BAD MOJO: exceeded memory requirement in extract_filter **\n");
		exit(-1);
	      }
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


// *************************** GAUSSIAN FILTERS ************************

void extract_filter_gaussian(double* rhoR, double* rhoI, 
			     int N[3], double latt[9],
			     double& l0) 
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
  int h[3], half[3];
  int MEMORY = (Nh[0]+1)*(Nh[1]+1)*(Nh[2]+1); // hardcoded, I know...
  double* G_l = new double[MEMORY];
  double* rho_l = new double[MEMORY];
  int Nl = 0;

  double rho_min = rhoR[0] * RHO_MIN_SCALE;
  
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

	if (scaled_rho < rho_min) continue;
        
	for (h[0]=0; h[0]<=half[0]; ++(h[0]))
	  for (h[1]=0; h[1]<=half[1]; ++(h[1]))
	    for (h[2]=0; h[2]<=half[2]; ++(h[2])) {
	      double nc_h[3];
	      for (int d=0; d<3; ++d) nc_h[d] = nc[d] - h[d]*N[d];
	      double Gabs2 = magnsq(Gmetric, nc_h);
	      add_list(G_l, rho_l, Nl, Gabs2, log(scaled_rho));
	      if (Nl >= MEMORY) {
		fprintf(stderr, "** BAD MOJO: exceeded memory requirement in extract_filter **\n");
		exit(-1);
	      }
	    }
      }

  // now, fit to a Gaussian (which is a linear fit, in this case)
  double G_sum=0, G2_sum=0, rho_sum=0, rhoG_sum=0;
  for (int i=0; i<Nl; ++i) {
    G_sum    += G_l[i];
    G2_sum   += G_l[i] * G_l[i];
    rho_sum  += rho_l[i];
    rhoG_sum += rho_l[i] * G_l[i];
  }
  // garbage collection
  delete[] G_l;
  delete[] rho_l;

  l0 = (Nl*rhoG_sum - G_sum*rho_sum)/(Nl*G2_sum - G_sum*G_sum);
  if (l0>0) {
    fprintf(stderr, "Somehow did not get a _decreasing_ Gaussian in extract_filter_gaussian\n");
    l0 = 0;
  }
  l0 = sqrt(-l0)*(2.*M_PI);  // scale accordingly.
}


// assumes that l0 has been scaled by some factor
void apply_filter_gaussian(double* rhoR, double* rhoI, 
			   int N[3], double latt[9],
			   double l0) 
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
  double lsquare = (l0*l0)/(4.*M_PI*M_PI);
  
  for (n[2]=0; n[2]<N[2]; ++(n[2]))
    for (n[1]=0; n[1]<N[1]; ++(n[1]))
      for (n[0]=0; n[0]<N[0]; ++(n[0])) {
	for (int d=0; d<3; ++d) {
	  if (n[d]>Nh[d]) nc[d] = n[d]-N[d];
	  else            nc[d] = n[d];
	}

	double Gabs2 = magnsq(Gmetric, nc);
	double scale = exp(-Gabs2*lsquare);
	rhoR[nind] *= scale;
	rhoI[nind] *= scale;
	++nind;
      }
}
