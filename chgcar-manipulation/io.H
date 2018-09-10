#ifndef __IO_H
#define __IO_H

/*
  Program: io.H
  Author:  D. Trinkle
  Date:    Jan. 22, 2003
  Purpose: This is the general part of the code, designed to do all
           parsing of both the command line and the (ohmms style)
	   metafiles.

	   This header file contains the following subroutines:

	   parse_commandline (argc, argv, Narg, args, 
                              VERBOSE, TESTING, MEMORY,
			      Nflags, flaglist[], flagon[]);
	     --separates the commandline arguments given in argv
	       into the Narg pieces (in order), and removes any
	       verbose, testing, or memory flags and puts them into the
	       associated flags.  We also check flags against a user
	       defined list of possible flags in flaglist, and
	       turn on the appropriate flagon if found.  Doesn't
	       support flags with parameters.  Return value: ERROR code.
	       Note: there is a special error code, for help (ERROR_HELP).

	   read_meta (meta_file, cart, u, Natoms)

	     --reads in the input meta file (formatted for OHMMS).
               Assumes that meta_file is an already open file (so that
               we can feed stdin as an option).  Doesn't close the
               file either.  Allocates memory for each of the
               pointers, and fills it up.  Return value: ERROR code.

           free_meta (u, Natoms)
	     --frees up the associated memory (done for completeness only).

	   verbose_output_meta (cart, u, Natoms, opening) 
	     --outputs everything we read in from the input file,
	       preceding each line with opening.

	   write_meta (meta_file, cart, u, atom_type, Natoms)
	     --creates a meta_file with the correctly labeled atoms

	   write_xyz (meta_file, cart, u, atom_type, Natoms, basename,
	              Nexpand, expand)
	     --creates an xyz file with the correctly labeled atoms.
	       We name each atom as "basename" with the atom_type
	       number concatenated.  Nexpand is the number of "expansion"
	       integer shifts we give (we output Natoms * Nexpand atoms).


	   print_long_help (prog_name, arglist, Nflags, flaglist, argexpl)
	     --prints a standard help message (includes input file
	       format, list of crystal classes), but also lists
	       the correct format for the program by calling
	       print_short_help().

	   print_short_help (prog_name, arglist, Nflags, flaglist, argexpl)
	     --lists the correct format for the program.  Called if
	       the number of arguments is too low, or by long_help.

	   insidecell (double x)
	     --translate x to within [0,1)

	   myopenr (char *filename)
	     --open a file; if filename is '-', set to stdin

	   myopenw (char *filename)
	     --open a file; if filename is '-', set to stdout

	   myclose (FILE *f)
	     --close a file, unless its stdin or stdout
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <iomanip>
#include "dcomp.H"

//****************************** ERROR FLAGS ***************************

// Note: we reuse error codes, since once we throw an error, everything
//       will stop anyway.

// errors from parse_commandline:
const int ERROR_SHORT   = 1; // Not enough arguments on line
const int ERROR_BADFLAG = 2; // We got a bad flag
const int ERROR_HELP    = 4; // Not really an error, but we can't continue

// errors from read_meta:
const int ERROR_NOFILE =      1; // General error if a file can't be opened.
const int ERROR_BADFILE =     2; // File has wrong format.
const int ERROR_ZEROVOL =     4; // Cartesian cell has zero volume
const int ERROR_LEFTHANDED =  8; // Cartesian cell has negative volume
const int ERROR_MEMORY =     16; // Not enough memory


inline int has_error (int ERROR, int ERROR_CODE) 
{
  return ( (ERROR & ERROR_CODE) != 0 );
}


//****************************** SUBROUTINES ***************************
// Headers for all:
inline double insidecell (double x);
int parse_commandline (int argc, char **argv, 
		       int Narg, char **args, 
		       int &VERBOSE, int &TESTING, int &MEMORY,
		       int Nflags, const char* flaglist, int* flagon);

// handle a variable number of arguments... (called by other)
int parse_commandline_var (int argc, char **argv, 
			   int& Nargs, char **args, 
			   int &VERBOSE, int &TESTING, int &MEMORY,
			   int Nflags, const char* flaglist, int* flagon);

int  read_meta (FILE *meta_file, double cart[9], double** &u, int &Natoms);
void free_meta (double** &u, int Natoms);
void verbose_output_meta (double cart[9], double** u, int Natoms,
			  const char* opening);
void write_meta (FILE *meta_file, double cart[9], double** u,
		 int* atom_type, int Natoms, char* basename, 
		 char** phase_name, int output_mask);
void write_xyz (FILE *meta_file, double cart[9], double** u,
		int* atom_type, int Natoms, char* basename,
		char** phase_name, int output_mask, 
		int Nexpand, int** expand);

void print_short_help (char *prog_name, const char *arglist, 
		       int Nflags, const char* flaglist, const char *argexpl);
void print_long_help (char *prog_name, const char *arglist,
		      int Nflags, const char* flaglist, const char *argexpl);


inline FILE* myopenr(char *filename) 
{
  if (filename[0] == '-') return stdin;
  else                    return fopen(filename, "r");
}

inline FILE* myopenw(char *filename) 
{
  if (filename[0] == '-') return stdout;
  else                    return fopen(filename, "w");
}
  
inline void myclose(FILE *f) 
{
  if ((f != stdin) && (f != stdout) && (f != NULL)) fclose(f);
}


//******************************* insidecell ***************************
inline double insidecell (double x) 
{
  double y;
  // Return a value less than 1, and greater than or equal to 0.
  for (y = x + 1 + 1.e-13 - (int)(x); y >= 1.0; y -= 1.0) ;
  y -= 1.e-13;
  if (y < 0.) return 0;
  return y;
}


// ***************************** is_blank ******************************
inline int is_blank (char c) {return (c == ' ') || (c == '\t') || (c == '\n');}


// Move to the next non-blank part of string
inline void nextnonblank (char* &p) 
{
  for ( ; ((*p)!='\0') && is_blank(*p); ++p) ;
}

// Move to the end of string (find the first blank)
inline void nextblank (char* &p) 
{
  for ( ; ((*p)!='\0') && (!is_blank(*p)); ++p) ;
}

const char COMMENT_CHAR = '#';

inline void nextnoncomment (FILE* infile, char* dump, const int &size) 
{
  do {fgets(dump, size, infile);
  } while ((!feof(infile)) && (dump[0] == COMMENT_CHAR));
}


//*************************** parse_commandline ************************
// We separate the commandline arguments given in argv into the Narg 
// pieces (in order), and remove verbose, testing, or memory flags and
// puts them into the associated flags.  Return value: ERROR code.
// Note the error codes listed above.
// 
// The flags (not case sensitive) -- anything starting with -, but not
// a dash by itself (could refer to stdin / stdout)
//   -v       : verbose
//   -t       : testing mode (also flips verbose, by default)
//   -m <mem> : a general memory flag; set MEMORY to value of second param.
//   -h       : we _want_ help
//   (all other flags generate an error)
// Everything that's _not_ a flag gets put into args, one at a time.
// We assume that args has NOT BEEN ALLOCATED, and that
// VERBOSE, TESTING, and MEMORY have already been set to default
// values--we will only change them if they show up in the command line.
// 
// We are going to do pointer _assignments_ on args, so it should
// simply be an array of N pointers to char* (there's no reason to 
// go around _copying_ these arguments).

// Now backwards compatible; the _var version allows for a variable
// number of arguments (-1), but still throws errors like a champ
// if given a constant number of arguments

int parse_commandline(int argc, char **argv, 
		      int Narg, char **args, 
		      int &VERBOSE, int &TESTING, int &MEMORY,
		      int Nflags, const char* flaglist, int* flagon) 
{
  int N, ERROR;
  N=Narg;
  ERROR = parse_commandline_var(argc, argv, N, args, VERBOSE,
				TESTING, MEMORY, Nflags, flaglist, flagon);
  if (N<Narg)
    ERROR = ERROR | ERROR_SHORT;
  return ERROR;
}

int parse_commandline_var(int argc, char **argv, 
			  int& Narg, char **args, 
			  int &VERBOSE, int &TESTING, int &MEMORY,
			  int Nflags, const char* flaglist, int* flagon) 
{
  int ERROR, HELP;
  int count; // Count of how many args we have.
  int n;     // Counting through argc.
  int i, j;  // Counting through our flaglist.
  char flag; // Our flag (maybe)
  int FLAGFOUND;
  int MEMFOUND; // Special flag indicating that we hit -m
  
  count = 0; // No arguments yet.
  ERROR = 0; // We just started, so nothing's gone wrong yet.  
  HELP = 0;  // Assume they don't want help till they ask for it.
  // Note: we DON'T set verbose or testing... if they're already
  // true before we get to this point, we don't need to upset that.

  // First argument is always the program name, so start at 1.
  for (n=1; n<argc; ++n) {
    if ((argv[n][0] == '-') && (argv[n][1] != 0) && (!isdigit(argv[n][1]))) {
      // We've got a flag!!
      MEMFOUND = 0;
      for (j=1; argv[n][j] != 0; ++j) {
	flag = tolower(argv[n][j]);
	switch (flag) {
	case 'v': VERBOSE = 1; break;
	case 't': TESTING = 1; break;
	case 'h': HELP = 1; break;
	case 'm': MEMFOUND = 1; break; // We'll extract MEMORY later...
	  // If we get this far, that means we may have a user-defined
	  // flag
	default:
	  FLAGFOUND = 0;
	  for (i=0; i<Nflags; ++i)
	    if (flag == tolower(flaglist[i])) {
	      flagon[i] = 1;
	      FLAGFOUND = -1;
	    }
	  if (!FLAGFOUND)
	    // Throw an error--we didn't find our flag.
	    ERROR = ERROR | ERROR_BADFLAG;
	}
      }
      if (MEMFOUND) {
	// Extract the "memory" value:
	if ((n+1)<argc) {
	  ++n; // Increment arg counter:
	  sscanf(argv[n], "%d", &MEMORY);
	}
	else
	  ERROR = ERROR | ERROR_SHORT;
      }
    }
    else {
      // We didn't get a flag, so copy over the arg. info, if 
      // we haven't read all the arguments yet.
      // Narg < 0 --> variable number of arguments
      if ( (Narg<0) || (count < Narg) ) {
	// POINTER ASSIGNMENT:
	args[count] = argv[n];
	++count;
      }
    }
  }
  // If we're testing, we're _definitely_ verbose:
  VERBOSE = VERBOSE || TESTING;
  if (HELP)         ERROR = ERROR | ERROR_HELP;  // We called help.
  if (count < abs(Narg)) ERROR = ERROR | ERROR_SHORT; // Not enough args.
  Narg = count;
  
  return ERROR;
}


//****************************** read_meta *****************************
// We read in the input file (formatted as an OHMMS metafile).  Assumes
// that meta_file is an already open file (so that we can feed stdin as
// an option).  Doesn't close the file either.
// Note: an OHMMS meta_file contains lots of information; we're looking
// for just two tags: <lattice>, for the lattice information, and
// <TBMD>, for the atom location information.  We ignore everything else.
// Also assumes the atoms are in unit coord. (though this could be
// changed for future versions to convert back to unit cell coord).
// Allocates memory for each of the pointers, and fills it up.
// Return value: ERROR code (see above)
// ==== input file ====
// ...                        # extra stuff (we ignore)
// <lattice>
// a1.x a1.y a1.z             # Cartesian coord of unit cell
// a2.x a2.y a2.z
// a3.x a3.y a3.z
// </lattice>
// ...                        # extra stuff
// <TBMD>
// number_of_ions = Natoms    # Number of atoms in unit cell
// positions unit/cartesian   # Default is cartesian, unless "u" is first
// i1 u1.1 u1.2 u1.3          # Atom type, and locations in direct coord.
// ...
// iN uN.1 uN.2 uN.3
// </TBMD>
// ...                        # extra stuff
// ==== input file ====

int read_meta(FILE *meta_file, double cart[9], double** &u, int &Natoms) 
{
  char dump[512]; // For reading files (allows extra characters at the end
                  // of each line)
  int i, k;
  int READ_LATTICE, READ_ATOMS;
  int IN_UNIT;
  char *p;        // String used to "parse" our tags

  if (meta_file == NULL)
    return ERROR_NOFILE;

  // ========================== read input file ========================

  // We're gonna read until either:
  // a) we've read both our lattice and atoms, or
  // b) EOF
  READ_LATTICE = 0;
  READ_ATOMS = 0;
  while ( ( (!READ_LATTICE) || (!READ_ATOMS) )
	  && (! feof(meta_file)) ) {
    fgets(dump, sizeof(dump), meta_file); // Grab our next line.
    // First, convert dump to lower case:
    for (p = dump; *p != 0; ++p) *p = tolower(*p);

    // <lattice>
    p = strstr(dump, "lattice");
    if ((p != NULL) && (!READ_LATTICE)) {
      // We've got lattice info!
      READ_LATTICE = -1;
      // a1.x a1.y a1.z                # Cartesian coord of unit cell
      // a2.x a2.y a2.z
      // a3.x a3.y a3.z
      for (i=0; i<3; ++i) {
	fgets(dump, sizeof(dump), meta_file);
	sscanf(dump, "%lf %lf %lf", 
	       cart+(3*i), cart+(3*i+1), cart+(3*i+2));
      }
      continue; // Go on to the next line in the file.
    }

    // number_of_ions (easier to search for than TBMD)
    p = strstr(dump, "number_of_ions");
    if ((p != NULL) && (!READ_ATOMS)) {
      // We've got atom info!
      READ_ATOMS = -1;

      // number_of_ions = Natoms
      // Let's jump to the first integer on that line
      p = dump + (strcspn (dump, "0123456789"));
      sscanf(p, "%d", &Natoms);

      // positions unit/cartesian  # Read that second word.
      fgets(dump, sizeof(dump), meta_file);
      p = dump; nextnonblank(p); // First word
      // Now, let's look for a word beginning with u (after "positions"):
      do {
	nextblank(p); nextnonblank(p); // Next word
	IN_UNIT = ((*p == 'u') || (*p == 'U'));
      } while ((!IN_UNIT) && (*p != '\0') && (*p != '#'));
      
      // i1 u1.1 u1.2 u1.3          # Atom type, and locations in direct coord.
      // ...
      // iN uN.1 uN.2 uN.3
      u = new double *[Natoms];
      for (i=0; i<Natoms; ++i) {
	u[i] = new double[3];
	fgets(dump, sizeof(dump), meta_file);
	// Note: we ignore the "atom type" in the metafile.
	sscanf(dump, "%d %lf %lf %lf", &k, u[i], u[i]+1, u[i]+2);
	if (IN_UNIT) {
	  for (k=0; k<3; ++k)
	    u[i][k] = insidecell(u[i][k]);
	}
	else {
	  // Do nothing... we'll do the conversion *after*, in case
	  // we haven't read our cell information yet.
	}
      }
      continue;
    }
    // Space reserved for possible future parsing
  }

  if ( (!READ_LATTICE) || (!READ_ATOMS) ) // We didn't get everything
    return ERROR_BADFILE;

  // ============================ test values ==========================

  // Determine determinants:
  double det_cart;
  
  det_cart = cart[0]*(cart[4]*cart[8] - cart[7]*cart[5]) -
    cart[1]*(cart[3]*cart[8] - cart[6]*cart[5]) +
    cart[2]*(cart[3]*cart[7] - cart[6]*cart[4]);
  
  if (dcomp(det_cart, 0)) return ERROR_ZEROVOL;
  if (det_cart < 0)       return ERROR_LEFTHANDED;

  // Now, convert cartesian to unit cell if necessary:
  if (! IN_UNIT) {
    // First, invert the cell.
    double inv_cart[9];
    double unit[3];
    inv_cart[0] = cart[4]*cart[8] - cart[5]*cart[7];
    inv_cart[1] = cart[2]*cart[7] - cart[1]*cart[8];
    inv_cart[2] = cart[1]*cart[5] - cart[2]*cart[4];
    
    inv_cart[3] = cart[5]*cart[6] - cart[3]*cart[8];
    inv_cart[4] = cart[0]*cart[8] - cart[2]*cart[6];
    inv_cart[5] = cart[2]*cart[3] - cart[0]*cart[5];
    
    inv_cart[6] = cart[3]*cart[7] - cart[4]*cart[6];
    inv_cart[7] = cart[1]*cart[6] - cart[0]*cart[7];
    inv_cart[8] = cart[0]*cart[4] - cart[1]*cart[3];
    for (i=0; i<9; ++i) inv_cart[i] *= (1./det_cart);

    // Now, let's convert 'em:
    for (i=0; i<Natoms; ++i) {
      for (k=0; k<3; ++k)
	unit[k] = insidecell(u[i][0]*inv_cart[k] +
			     u[i][1]*inv_cart[3+k] +
			     u[i][2]*inv_cart[6+k]);
      // Now, copy over.
      for (k=0; k<3; ++k) u[i][k] = unit[k];
    }
  }

  // If we've gotten this far, it all went well.
  return 0;
}


// ************************* verbose_output_meta ***********************
// Little routine to output everything we read in from cell.

void verbose_output_meta (double cart[9], double** u, int Natoms,
			  const char* opening) 
{
  int i;

  printf("%sCartesian coordinates:\n", opening);
  for (i=0; i<3; ++i)
    printf("%sa%1d = %8.5lf %8.5lf %8.5lf\n", opening, i+1, 
	   cart[3*i], cart[3*i+1], cart[3*i+2]);

  printf("%sAtoms in cell (%d):\n", opening, Natoms);
  if (u == NULL)
    printf("%satom positions not read.\n", opening);
  else
    for (i=0; i<Natoms; ++i)
      printf("%su%1d = %8.5lf %8.5lf %8.5lf\n", opening, i+1, 
	     u[i][0], u[i][1], u[i][2]);

}


//****************************** write_meta ****************************
// Output the metafile with the correct atom labels
void write_meta (FILE *meta_file, double cart[9], double** u,
		 int* atom_type, int Natoms, char* basename, 
		 char** phase_name, int output_mask) 
{
  int n, j;

  int omask, WITH_NONE;
  if (output_mask <= 0) {
    WITH_NONE = -1;
    omask = -output_mask;
  }
  else {
    WITH_NONE = 0;
    omask = output_mask;
  }
    
  // How many atoms?
  int count;
  int *show_atom;
  show_atom = new int[Natoms];

  for (n=0; n<Natoms; ++n) {
    show_atom[n] = 0;
    if (atom_type[n] == 0) {
      if (WITH_NONE) show_atom[n] = 1;
    }
    else {
      for (j=0; j<30; ++j) {
	if ( ( (1<<j) & omask ) != 0)
	  if ( ( (1<<j) & atom_type[n]) != 0) show_atom[n] = 1;
      }
    }
  }
  for (n=0, count=0; n<Natoms; ++n)
    count += show_atom[n];

  fprintf(meta_file, "<unitcell>\n");
  fprintf(meta_file, "  <lattice>\n");
  fprintf(meta_file, "    %19.15lf %19.15lf %19.15lf  # a1\n",
	  cart[0], cart[1], cart[2]);
  fprintf(meta_file, "    %19.15lf %19.15lf %19.15lf  # a2\n",
	  cart[3], cart[4], cart[5]);
  fprintf(meta_file, "    %19.15lf %19.15lf %19.15lf  # a3\n",
	  cart[6], cart[7], cart[8]);
  fprintf(meta_file, "  </lattice>\n");
  fprintf(meta_file, "  <bconds> p p p </bconds>\n");
  fprintf(meta_file, "</unitcell>\n");
  fprintf(meta_file, "<ptclset>\n");
  fprintf(meta_file, "  <format> tbmd </format>\n");
  fprintf(meta_file, "  <name> none </name>\n");
  fprintf(meta_file, "<TBMD>\n");
  fprintf(meta_file, "number_of_ions = %d\n", count);
  fprintf(meta_file, "positions = unit\n");

  for (n=0; n<Natoms; ++n)
    if (show_atom[n]) {
      fprintf(meta_file, "  1  %17.15lf %17.15lf %17.15lf  # %d %s",
	      u[n][0], u[n][1], u[n][2], n+1, basename);
      for (j=0; j<30; ++j)
	if ( ( (1<<j) & atom_type[n]) != 0)
	  fprintf(meta_file, " %s", phase_name[j]);
      fprintf(meta_file, "\n");
    }
  
  fprintf(meta_file, "</TBMD>\n");
  fprintf(meta_file, "</ptclset>\n");
}


//****************************** write_xyz *****************************
// Output an XYZ format file.  We make "copies" of the atoms by using
// the "expand" vectors.  If Nexpand == 0, we just output what we've
// got; else, we'll output Natoms * Nexpand atoms.
void write_xyz (FILE *meta_file, double cart[9], double** u,
		int* atom_type, int Natoms, char* basename,
		char** phase_name, int output_mask, 
		int Nexpand, int** expvect)
{
  int n, i, j;
  double xyz[3];
  double* uvect;
  int expand_out;
  int nexpand, **expand;
  
  expand_out = ((Nexpand > 1) && (expvect != NULL));
  if (expand_out) {
    nexpand = Nexpand;
    expand = expvect;
  }
  else {
    nexpand = 1;
    expand = new int *[1];
    expand[0] = new int[3];
    for (i=0; i<3; ++i) expand[0][i] = 0;
  }

  int omask, WITH_NONE;
  if (output_mask <= 0) {
    WITH_NONE = -1;
    omask = -output_mask;
  }
  else {
    WITH_NONE = 0;
    omask = output_mask;
  }
    
  // How many atoms?
  int count;
  int *show_atom;
  show_atom = new int[Natoms];

  for (n=0; n<Natoms; ++n) {
    show_atom[n] = 0;
    if (atom_type[n] == 0) {
      if (WITH_NONE) show_atom[n] = 1;
    }
    else {
      for (j=0; j<30; ++j) {
	if ( ( (1<<j) & omask ) != 0)
	  if ( ( (1<<j) & atom_type[n]) != 0)
	    show_atom[n] = 1;
      }
    }
  }
  for (n=0, count=0; n<Natoms; ++n)
    count += show_atom[n];

  fprintf(meta_file, "%d\n", count*nexpand);
  fprintf(meta_file, "  %s phases: ", basename);
  for (j=0; j<30; ++j)
    if ( ( (1<<j) & output_mask ) != 0)
      fprintf(meta_file, "%s ", phase_name[j]);
  if (WITH_NONE)
    fprintf(meta_file, "unidentified ");
  fprintf(meta_file, ": flagged, expanded %d times\n", nexpand);

  for (n=0; n<Natoms; ++n)
    if (show_atom[n]) {
      uvect = u[n];
      for (i=0; i<nexpand; ++i) {
	xyz[0] = cart[0]*(uvect[0] + expand[i][0])
	  + cart[3]*(uvect[1] + expand[i][1])
	  + cart[6]*(uvect[2] + expand[i][2]);
	xyz[1] = cart[1]*(uvect[0] + expand[i][0])
	  + cart[4]*(uvect[1] + expand[i][1])
	  + cart[7]*(uvect[2] + expand[i][2]);
	xyz[2] = cart[2]*(uvect[0] + expand[i][0])
	  + cart[5]*(uvect[1] + expand[i][1])
	  + cart[8]*(uvect[2] + expand[i][2]);
	fprintf(meta_file, "%s", basename);
	for (j=0; j<30; ++j)
	  if ( ( (1<<j) & atom_type[n]) != 0)
	    fprintf(meta_file, "_%s", phase_name[j]);
	fprintf(meta_file, " %19.15lf %19.15lf %19.15lf\n",
		xyz[0], xyz[1], xyz[2]);
      }
    }
  
  if (!expand_out) {
    delete[] expand[0];
    delete[] expand;
  }
  else expand = NULL;

  delete[] show_atom;
}



//****************************** free_meta *****************************
// Free up the group and atom info:

void free_meta(double** &u, int Natoms) 
{
  int i;

  if (u != NULL) {
    for (i=0; i<Natoms; ++i)
      delete[] u[i];
    delete[] u;
  }
}



// ************************** print_short_help *************************
// List the correct format for the program.  Called if the number of
// arguments is too low, or by print_long_help.
// Prog_name should argv[0] (the invoked program name)
// arglist does NOT need to contain the flag info
// argexpl can contain an optional explanation of the command line param

const char* FLAGLIST = "[-v] [-t] [-m <memorysize>] [-h]";
const char* FLAGexpl = 
"  -v               verbose output\n\
  -t               testing (sets -v as well)\n\
  -m <memorysize>  set memory allocation to <memorysize>\n\
  -h               print detailed help";

void print_short_help (char *prog_name, const char *arglist, 
		       int Nflags, const char* flaglist, const char *argexpl)
{
  int i;
  fprintf(stderr, "%s %s", prog_name, arglist);
  if (Nflags > 0) fprintf(stderr, " [-");
  for (i=0; i<Nflags; ++i)
    fprintf(stderr, "%c", flaglist[i]);
  if (Nflags > 0) fprintf(stderr, "]");
  fprintf(stderr, " %s\n", FLAGLIST);
  if (argexpl != NULL)
    fprintf(stderr, "%s\n", argexpl);
  fprintf(stderr, "%s\n", FLAGexpl);
  fprintf(stderr, "\n");
}



// *************************** print_long_help *************************
// Print a standard help message (includes input file format),
// but also lists the correct format for the program by calling
// print_short_help().

const char* INPUTFILE = 
"  We ignore all of the tags in the meta_file except:\n\
   <lattice>\n\
   a1.x a1.y a1.z             # Cartesian coord of unit cell\n\
   a2.x a2.y a2.z\n\
   a3.x a3.y a3.z\n\
\n\
   and:\n\
   <TBMD>\n\
   number_of_ions = Natoms    # Number of atoms in unit cell\n\
   positions unit/cartesian   # Default is cartesian, unless \"u\" is first\n\
   i1 u1.1 u1.2 u1.3          # Atom type, and locations in direct coord.\n\
   ...\n\
   iN uN.1 uN.2 uN.3";

void print_long_help (char *prog_name, const char *arglist,
		      int Nflags, const char* flaglist, const char *argexpl)
{
  // Call short help first.
  print_short_help(prog_name, arglist, Nflags, flaglist, argexpl);
  fprintf(stderr, "OHMMS meta-file format:\n%s\n\n", INPUTFILE);
}
#endif
