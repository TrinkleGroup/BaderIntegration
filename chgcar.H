#ifndef __CHGCAR_H
#define __CHGCAR_H
#include <math.h>
#include <string.h>

/*
  Program: chgcar.H
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: Procedures for dealing with reading chgcar files.
*/


void fort_printf(char* str, const double& x);

void read_grid(FILE* infile, double* rho, int Ng);
void read_grid(FILE* infile, double* rho, int Ng, int* index);

void write_grid(FILE* outfile, double* rho, int Ng, int* index);
void write_grid(FILE* outfile, double* rho, int Ng);

void read_CHGCAR_header(FILE* infile, char* comment, double& a0,
			double latt[9], int* &Natom, double** &uatom,
			int READATOMS);

void write_CHGCAR_header(FILE* outfile, char* comment, double a0,
			 double latt[9], int* Natom, double** uatom);

void skip_CHGCAR_header(FILE* infile);


/*================================= main ==================================*/

// this is how many entries appear per line...
const int READ_stride = 5;


void read_grid(FILE* infile, double* rho, int Ng) 
{
  if (infile==NULL) return;

  int n, line, nl;
  char dump[512];

  if (Ng == 0) {
    fgets(dump, sizeof(dump), infile);
    return;
  }

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
    char *strp, *endp;
    for (strp=dump; n<Ng; ++n, strp=endp)
      rho[n] = strtod(strp, &endp);
    strp=NULL;
    endp=NULL;
  }
}

// allows an "index" function to be specified for sorting
void read_grid(FILE* infile, double* rho, int Ng, int* index) 
{
  if (infile==NULL) return;

  int n, line, nl;
  char dump[512];

  if (Ng == 0) {
    fgets(dump, sizeof(dump), infile);
    return;
  }

  nl = Ng/READ_stride;
  n=0;
  
  for (line=0, n=0; line<nl; ++line, n+=READ_stride) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf %lf %lf",
	   rho+index[n], rho+index[n+1], rho+index[n+2], rho+index[n+3], rho+index[n+4]);
  }

  if (n!=Ng) {
    // last line...
    fgets(dump, sizeof(dump), infile);
    char *strp, *endp;
    for (strp=dump; n<Ng; ++n, strp=endp)
      rho[index[n]] = strtod(strp, &endp);
    strp=NULL;
    endp=NULL;
  }
}

// outputs in true fortran style... we seem to be having trouble
// because we weren't using this... *sigh*
void fort_printf(char* str, const double& x) 
{
  sprintf(str, "%+17.10lE", x);
  str[2] = str[1];
  str[1] = '.';
  if (x >= 0) str[0] = '0';
  // fix the exponent (only if non-zero)
  if (x != 0) {
    int expon = strtol(str+14, (char**)NULL, 10);
    expon++;
    if (expon>=0) str[14] = '+';
    else         {str[14] = '-'; expon = -expon;}
    sprintf(str+15, "%02d", expon);
  }
}



// not quite the same as what is output by CHGCAR, but maintains the same
// accuracy, so "just as good", I'd say.
void write_grid(FILE* outfile, double* rho, int Ng) 
{
  if (outfile==NULL) return;

  if (Ng == 0) {
    fprintf(outfile, "\n");
    return;
  }
  int n;
  char dump[18];
  for (n=0; n<Ng; ++n) {
    fort_printf(dump, rho[n]);
    fprintf(outfile, " %s", dump);
    if (n%5 == 4) fprintf(outfile, "\n");
  }
  if (Ng%5 != 0) fprintf(outfile, "\n");
}

// allows an "index" function to be specified for sorting
void write_grid(FILE* outfile, double* rho, int Ng, int* index)
{
  if (outfile==NULL) return;

  if (Ng == 0) {
    fprintf(outfile, "\n");
    return;
  }
  int n;
  char dump[18];
  for (n=0; n<Ng; ++n) {
    fort_printf(dump, rho[index[n]]);
    fprintf(outfile, " %s", dump);
    if (n%5 == 4) fprintf(outfile, "\n");
  }
  if (Ng%5 != 0) fprintf(outfile, "\n");
}


inline void read_CHGCAR_header(FILE* infile, char* comment, double& a0,
			       double latt[9], int* &Natom, double** &uatom) 
{
  read_CHGCAR_header(infile, comment, a0, latt, Natom, uatom, -1);
}

inline void read_CHGCAR_header(FILE* infile, char* comment, double& a0,
			       double latt[9]) 
{
  int* Natom=NULL;
  double** uatom=NULL;
  read_CHGCAR_header(infile, comment, a0, latt, Natom, uatom, 0);
}

void read_CHGCAR_header(FILE* infile, char* comment, double& a0,
			double latt[9], int* &Natom, double** &uatom,
			int READATOMS) 
{
  if (infile==NULL) return;

  char dump[512];
  int i;

  //++ comment line
  fgets(dump, sizeof(dump), infile);
  char *p = 
    stpcpy(comment, dump);
  (p-1)[0] = '\0'; // remove the newline at the end of dump.

  //++ a0
  fgets(dump, sizeof(dump), infile);
  sscanf(dump, "%lf", &a0);

  //++ lattice vectors
  for (i=0; i<3; ++i) {
    fgets(dump, sizeof(dump), infile);
    sscanf(dump, "%lf %lf %lf", latt+i, latt+3+i, latt+6+i);
  }

  //++ number of atoms
  fgets(dump, sizeof(dump), infile);
  char *strp, *endp=dump;
  int types = 0, Natoms=0;
  if (READATOMS && (Natom==NULL)) Natom=new int[128];
  do {
    strp=endp;
    i = strtol(strp, &endp, 10);
    if (strp!=endp) {
      Natoms += i;
      if (READATOMS) Natom[types] = i;
      ++types;
    }
  } while (strp!=endp);
  if (READATOMS) Natom[types] = -1;

  //++ direct coord. (CHGCAR should always be output with "Direct")
  fgets(dump, sizeof(dump), infile);

  //++ atom positions
  if (READATOMS) uatom = new double*[Natoms];
  for (i=0; i<Natoms; ++i) {
    fgets(dump, sizeof(dump), infile);
    if (READATOMS) {
      uatom[i] = new double[3];
      double* in_vect = uatom[i];
      sscanf(dump, "%lf %lf %lf", in_vect, in_vect+1, in_vect+2);
    }
  }
  
  //++ blank line
  fgets(dump, sizeof(dump), infile);
}


void skip_CHGCAR_header(FILE* infile) 
{
  if (infile==NULL) return;

  char dump[512];

  //++ comment line
  fgets(dump, sizeof(dump), infile);

  //++ a0
  fgets(dump, sizeof(dump), infile);

  //++ lattice vectors
  for (int i=0; i<3; ++i) fgets(dump, sizeof(dump), infile);

  //++ number of atoms
  fgets(dump, sizeof(dump), infile);
  char *strp, *endp=dump;
  int Natoms=0;
  do {
    strp=endp;
    int i = strtol(strp, &endp, 10);
    if (strp!=endp) Natoms += i;
  } while (strp!=endp);

  //++ direct coord. (CHGCAR should always be output with "Direct")
  fgets(dump, sizeof(dump), infile);

  //++ atom positions
  for (int i=0; i<Natoms; ++i) 
    fgets(dump, sizeof(dump), infile);
  
  //++ blank line
  fgets(dump, sizeof(dump), infile);
}


void write_CHGCAR_header(FILE* outfile, char* comment, double a0,
			 double latt[9], int* Natom, double* uatom[3]) 
{
  if (outfile==NULL) return;

  fprintf(outfile, "%s\n", comment);
  fprintf(outfile, "%20.17lf\n", a0);
  for (int i=0; i<3; ++i)
    fprintf(outfile, "  %11.6lf %11.6lf %11.6lf\n", latt[i], latt[i+3], latt[i+6]);
  int Natoms=0;
  for (int n=0; Natom[n]>=0; ++n) {
    Natoms += Natom[n];
    if (Natom[n]>0) fprintf(outfile, " %3d", Natom[n]);
  }
  fprintf(outfile, "\nDirect\n");
  for (int n=0; n<Natoms; ++n)
    fprintf(outfile, " %9.6lf %9.6lf %9.6lf\n", uatom[n][0],uatom[n][1],uatom[n][2]);
  fprintf(outfile, "\n");
}

#endif
