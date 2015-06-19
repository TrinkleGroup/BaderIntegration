#ifndef __MATRIX_H
#define __MATRIX_H
#include <math.h>

/*
  Program: matrix.H
  Author:  D. Trinkle
  Date:    Sept. 27, 2002
  Purpose: Define a series of matrix manipulations, all inlined, for
           3 x 3 matrices.  Designed to allow matrix manipulation
	   on integer and double matrices, using a 9-vector notation.

	   The nine vector index is defined as (i,j) -> (i*3 + j), or:

	   | 0 1 2 |
	   | 3 4 5 |
	   | 6 7 8 |

	   We first define some conversion utilities to take a matrix
	   M into it's 9-v notation, and also a quick indexing
	   function index(i,j) that produces the correct index, so
	   that v[index(i,j)] gives the correct member of the matrix.

	   We define: determinant, inverse (after dividing by the
	   determinant, which is returned), matrix multiplication,
	   matrix "squaring" (aT a), matrix comparison, rotation,
	   eigenvalues, and eigenvectors.
*/

const double ident[9] = {1,0,0,0,1,0,0,0,1};
const int ident_i[9] = {1,0,0,0,1,0,0,0,1};

// Convert between 3 x 3 matrix notation to 9 vector notation:
// We do all of our calculations, etc., in 9-v notation.
inline void MtoV (int m[3][3], int v[9]) 
{
  v[0] = m[0][0];  v[1] = m[0][1];  v[2] = m[0][2];
  v[3] = m[1][0];  v[4] = m[1][1];  v[5] = m[1][2];
  v[6] = m[2][0];  v[7] = m[2][1];  v[8] = m[2][2];
}

inline void VtoM (int v[9], int m[3][3])
{
  m[0][0] = v[0];  m[0][1] = v[1];  m[0][2] = v[2];
  m[1][0] = v[3];  m[1][1] = v[4];  m[1][2] = v[5];
  m[2][0] = v[6];  m[2][1] = v[7];  m[2][2] = v[8];
}

inline void MtoV (double m[3][3], double v[9]) 
{
  v[0] = m[0][0];  v[1] = m[0][1];  v[2] = m[0][2];
  v[3] = m[1][0];  v[4] = m[1][1];  v[5] = m[1][2];
  v[6] = m[2][0];  v[7] = m[2][1];  v[8] = m[2][2];
}

inline void VtoM (double v[9], double m[3][3])
{
  m[0][0] = v[0];  m[0][1] = v[1];  m[0][2] = v[2];
  m[1][0] = v[3];  m[1][1] = v[4];  m[1][2] = v[5];
  m[2][0] = v[6];  m[2][1] = v[7];  m[2][2] = v[8];
}

// Quick calculation to return the correct index in the 9-v
inline int index (int i, int j) 
{
  return 3*i + j;
}


// **** Matrix manipulation, all done in 9-v notation

// DET: Returns the det. of an integer matrix
inline int det (const int x[9]) 
{
  return (x[0]*(x[4]*x[8] - x[5]*x[7])
    + x[1]*(x[6]*x[5] - x[3]*x[8])
    + x[2]*(x[3]*x[7] - x[6]*x[4]));
}

inline double det (const double x[9]) 
{
  return (x[0]*(x[4]*x[8] - x[5]*x[7])
    + x[1]*(x[6]*x[5] - x[3]*x[8])
    + x[2]*(x[3]*x[7] - x[6]*x[4]));
}


// Rotate vectors in a matrix:
inline void rotate (const double a[9], double a_rot[9], int rot)
{
  int i;
  for (i=0; i<9; ++i)
    a_rot[i] = a[index(i/3,((i%3)+rot)%3)];
}

inline void rotate (const int a[9], int a_rot[9], int rot)
{
  int i;
  for (i=0; i<9; ++i)
    a_rot[i] = a[index(i/3,((i%3)+rot)%3)];
}


// Transpose of a matrix:
inline void transpose (const int a[9], int b[9])
{
  int i;
  for (i=0; i<9; ++i)
    b[i] = a[index(i%3, i/3)];
}

inline void transpose (const double a[9], double b[9])
{
  int i;
  for (i=0; i<9; ++i)
    b[i] = a[index(i%3, i/3)];
}

// Do the transpose "in place," with a loop "unrolled" version:

#define __ISWAP__(x,y) ({temp = x; x = y; y = temp;})
#define __DSWAP__(x,y) ({temp = x; x = y; y = temp;})

inline void self_transpose (int a[9]) 
{
  int temp;
  __ISWAP__(a[1], a[3]);
  __ISWAP__(a[2], a[6]);
  __ISWAP__(a[5], a[7]);
}

inline void self_transpose (double a[9]) 
{
  double temp;
  __DSWAP__(a[1], a[3]);
  __DSWAP__(a[2], a[6]);
  __DSWAP__(a[5], a[7]);
}

#undef __ISWAP__
#undef __DSWAP__


// MULT: evaluates the product of a * b, returns in c:
inline void mult (const int a[9], const int b[9], int c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

inline void mult (const double a[9], const int b[9], double c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

inline void mult (const int a[9], const double b[9], double c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

inline void mult (const double a[9], const double b[9], double c[9]) 
{
  c[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
  c[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
  c[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
  c[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
  c[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  c[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
  c[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
  c[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
  c[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

// MULT3 routines: D = ABC
// only dii is a little different... the rest are all the same
inline void mult (const int a[9], const int b[9], const int c[9],
		  int d[9]) 
{ int temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const int b[9], const int c[9],
		  double d[9]) 
{ int temp[9];  mult(b,c,temp);  mult(a, temp,d);}

inline void mult (const int a[9], const double b[9], const int c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const double b[9], const int c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const int a[9], const int b[9], const double c[9],
		  double d[9]) 
{ int temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const int b[9], const double c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const int a[9], const double b[9], const double c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}

inline void mult (const double a[9], const double b[9], const double c[9],
		  double d[9]) 
{ double temp[9];  mult(a,b,temp);  mult(temp,c,d);}



// MULT_VECT: Multiplies vector b by matrix a, result in c
inline void mult_vect (const int a[9], const int b[3], int c[3]) 
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
  c[2] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
}

inline void mult_vect (const double a[9], const int b[3], double c[3]) 
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
  c[2] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
}

inline void mult_vect (const int a[9], const double b[3], double c[3]) 
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
  c[2] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
}

inline void mult_vect (const double a[9], const double b[3], double c[3]) 
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
  c[2] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
}


// DOT PRODUCT:

inline int dot (const int a[3], const int b[3]) 
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double dot (const double a[3], const int b[3]) 
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double dot (const int a[3], const double b[3]) 
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double dot (const double a[3], const double b[3]) 
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


// MULT_VECT_MATRIX_VECT: Takes the inner product v1.a.v2

inline double innerprod (const int v1[3], const double a[9], const int v2[3])
{
  return ( v1[0]*a[0]*v2[0] + v1[0]*a[1]*v2[1] + v1[0]*a[2]*v2[2] +
           v1[1]*a[3]*v2[0] + v1[1]*a[4]*v2[1] + v1[1]*a[5]*v2[2] +
           v1[2]*a[6]*v2[0] + v1[2]*a[7]*v2[1] + v1[2]*a[8]*v2[2] );
}

inline double innerprod (const double v1[3], const double a[9], 
			 const int v2[3])
{
  return ( v1[0]*a[0]*v2[0] + v1[0]*a[1]*v2[1] + v1[0]*a[2]*v2[2] +
           v1[1]*a[3]*v2[0] + v1[1]*a[4]*v2[1] + v1[1]*a[5]*v2[2] +
           v1[2]*a[6]*v2[0] + v1[2]*a[7]*v2[1] + v1[2]*a[8]*v2[2] );
}

inline double innerprod (const int v1[3], const double a[9], 
			 const double v2[3]) 
{
  return ( v1[0]*a[0]*v2[0] + v1[0]*a[1]*v2[1] + v1[0]*a[2]*v2[2] +
	   v1[1]*a[3]*v2[0] + v1[1]*a[4]*v2[1] + v1[1]*a[5]*v2[2] +
	   v1[2]*a[6]*v2[0] + v1[2]*a[7]*v2[1] + v1[2]*a[8]*v2[2] );
}

inline double innerprod (const double v1[3], const double a[9], 
			 const double v2[3]) 
{
  return ( v1[0]*a[0]*v2[0] + v1[0]*a[1]*v2[1] + v1[0]*a[2]*v2[2] +
	   v1[1]*a[3]*v2[0] + v1[1]*a[4]*v2[1] + v1[1]*a[5]*v2[2] +
	   v1[2]*a[6]*v2[0] + v1[2]*a[7]*v2[1] + v1[2]*a[8]*v2[2] );
}


// CROSS_PRODUCT: evaluate a x b
inline void crossprod (const int a[3], const int b[3],
		       int c[3]) 
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

inline void crossprod (const double a[3], const int b[3],
		       double c[3]) 
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

inline void crossprod (const int a[3], const double b[3],
		       double c[3]) 
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

inline void crossprod (const double a[3], const double b[3],
		       double c[3]) 
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}



inline double tripleprod (const double a[3], const double b[3],
			  const double c[3]) 
{
  return c[0]*(a[1]*b[2] - a[2]*b[1])
    + c[1]*(a[2]*b[0] - a[0]*b[2])
    + c[2]*(a[0]*b[1] - a[1]*b[0]);
}


// MULT_SCALAR: Multiplies matrix a by scalar b, result in c
inline void mult (const int a[9], const int b, int c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const int b, const int a[9], int c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const int a[9], const double b, double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double b, const int a[9], double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double a[9], const int b, double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const int b, const double a[9], double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double a[9], const double b, double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }

inline void mult (const double b, const double a[9], double c[9]) 
{ for (int i=0; i<9; ++i) c[i] = b*a[i]; }


// SQUARE: Calculate xT * x
//   Note: the product is symmetric, so we get to do a smaller number
//   of multiplications.

inline void square (const int x[9], int s[9]) 
{
  s[0] = x[0]*x[0] + x[3]*x[3] + x[6]*x[6];
  s[1] = x[0]*x[1] + x[3]*x[4] + x[6]*x[7];
  s[2] = x[0]*x[2] + x[3]*x[5] + x[6]*x[8];
  s[3] = s[1];
  s[4] = x[1]*x[1] + x[4]*x[4] + x[7]*x[7];
  s[5] = x[1]*x[2] + x[4]*x[5] + x[7]*x[8];
  s[6] = s[2];
  s[7] = s[5];
  s[8] = x[2]*x[2] + x[5]*x[5] + x[8]*x[8];
}

inline void square (const int x[9], double s[9]) 
{
  s[0] = x[0]*x[0] + x[3]*x[3] + x[6]*x[6];
  s[1] = x[0]*x[1] + x[3]*x[4] + x[6]*x[7];
  s[2] = x[0]*x[2] + x[3]*x[5] + x[6]*x[8];
  s[3] = s[1];
  s[4] = x[1]*x[1] + x[4]*x[4] + x[7]*x[7];
  s[5] = x[1]*x[2] + x[4]*x[5] + x[7]*x[8];
  s[6] = s[2];
  s[7] = s[5];
  s[8] = x[2]*x[2] + x[5]*x[5] + x[8]*x[8];
}


inline void square (const double x[9], double s[9]) 
{
  s[0] = x[0]*x[0] + x[3]*x[3] + x[6]*x[6];
  s[1] = x[0]*x[1] + x[3]*x[4] + x[6]*x[7];
  s[2] = x[0]*x[2] + x[3]*x[5] + x[6]*x[8];
  s[3] = s[1];
  s[4] = x[1]*x[1] + x[4]*x[4] + x[7]*x[7];
  s[5] = x[1]*x[2] + x[4]*x[5] + x[7]*x[8];
  s[6] = s[2];
  s[7] = s[5];
  s[8] = x[2]*x[2] + x[5]*x[5] + x[8]*x[8];
}


// ********************************* magnsq ****************************
// Returns the squared magnitude of a vector, where metric[9] is the
// SYMMETRIC metric, and u is the vector in direct coordinates
inline double magnsq (const double metric[9], const double u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}

inline double magnsq (const int metric[9], const double u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}

inline double magnsq (const double metric[9], const int u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}

inline int magnsq (const int metric[9], const int u[3]) 
{
  return metric[0]*u[0]*u[0] + metric[4]*u[1]*u[1] + metric[8]*u[2]*u[2]
    + 2*metric[1]*u[0]*u[1] + 2*metric[2]*u[0]*u[2] + 2*metric[5]*u[1]*u[2];
}


// EQUAL: determine if two matrices are equal
inline int equal (const int a[9], const int b[9]) 
{
  return (a[0]==b[0]) && (a[1]==b[1]) &&(a[2]==b[2]) &&
    (a[3]==b[3]) && (a[4]==b[4]) && (a[5]==b[5]) &&
    (a[6]==b[6]) && (a[7]==b[7]) && (a[8]==b[8]);
}

// EQUAL_VECT: same thing for vectors
inline int equal_vect(const int a[3], const int b[3]) {
  return (a[0]==b[0]) && (a[1]==b[1]) && (a[2]==b[2]);
}

// ZERO_VECT
inline int zero_vect(const int a[3]) {
  return (a[0]==0) && (a[1]==0) && (a[2]==0);
}

// INVERSE: Evaluates inverse of x; returns it in inv / inverse (return value
// is det x).  Need to divide by det to get correct inverse.
inline int inverse (const int x[9], int inv[9])
{
  inv[0] = x[4]*x[8] - x[5]*x[7];
  inv[1] = x[2]*x[7] - x[1]*x[8];
  inv[2] = x[1]*x[5] - x[2]*x[4];

  inv[3] = x[5]*x[6] - x[3]*x[8];
  inv[4] = x[0]*x[8] - x[2]*x[6];
  inv[5] = x[2]*x[3] - x[0]*x[5];

  inv[6] = x[3]*x[7] - x[4]*x[6];
  inv[7] = x[1]*x[6] - x[0]*x[7];
  inv[8] = x[0]*x[4] - x[1]*x[3];
  
  return det(x);
}

inline double inverse (const double x[9], double inv[9])
{
  inv[0] = x[4]*x[8] - x[5]*x[7];
  inv[1] = x[2]*x[7] - x[1]*x[8];
  inv[2] = x[1]*x[5] - x[2]*x[4];

  inv[3] = x[5]*x[6] - x[3]*x[8];
  inv[4] = x[0]*x[8] - x[2]*x[6];
  inv[5] = x[2]*x[3] - x[0]*x[5];

  inv[6] = x[3]*x[7] - x[4]*x[6];
  inv[7] = x[1]*x[6] - x[0]*x[7];
  inv[8] = x[0]*x[4] - x[1]*x[3];
  
  return det(x);
}


// SYMM_INVERSE: inverts a symmetric matrix, keeping it symmetric
inline double symm_inverse (const double x[9], double inv[9])
{
  // enforce symmetrization:
  double x1 = 0.5*(x[1]+x[3]);
  double x2 = 0.5*(x[2]+x[6]);
  double x5 = 0.5*(x[5]+x[7]);
  double x0 = x[0];
  double x4 = x[4];
  double x8 = x[8];
  
  inv[0] = x4*x8 - x5*x5;
  inv[1] = x2*x5 - x1*x8;
  inv[2] = x1*x5 - x2*x4;

  inv[3] = inv[1];
  inv[4] = x0*x8 - x2*x2;
  inv[5] = x1*x2 - x0*x5;

  inv[6] = inv[2];
  inv[7] = inv[5];
  inv[8] = x0*x4 - x1*x1;

  // We roll our own slightly more stable determinant...
  // x0*x4*x8 - (x0*x5*x5 + x4*x2*x2 + x8*x1*x1)
  //  + 2.*x1*x2*x5;

  return 1./3.*(x0*inv[0]+x1*inv[1]+x2*inv[2] +
		x1*inv[1]+x4*inv[4]+x5*inv[5] +
		x2*inv[2]+x5*inv[5]+x8*inv[8]);
}


// CAREFUL_INVERSE: inverts a matrix, symmetric or not, but polishes
// the final result to within a set tolerance.
// NOTE: unlike all of our other inverses, we go ahead and scale by
// the determinant.
const double CAREFUL_INV_TOL = 1e-15;
const int MAX_ITER = 20; // maximum number of iterations

void careful_inverse (const double a[9], double inv[9], double &err) 
{
  double b[9], ba[9], bab[9];
  double deta;
  int d;
  
  deta = inverse(a, b); // our first crack at it...
  if (deta == 0) {
    err = 1; // a pretty obvious flag...
    return;
  }
  mult(b, 1./deta, b); // now b should be scaled accordingly.
  
  mult(b, a, ba);
  for (err=0, d=0; d<9; ++d) err += fabs(ba[d]-ident[d]);
  // polish the result, until err is below our desired tolerance
  // but don't do more than n steps.
  for (int n=0; (n<MAX_ITER) && (err > CAREFUL_INV_TOL); ++n) {
    mult(ba, b, bab);
    // polish in place:
    for (d=0; d<9; ++d) b[d] = 2.*b[d] - bab[d];
    mult(b, a, ba);
    for (err=0, d=0; d<9; ++d) err += fabs(ba[d]-ident[d]);
  }
  // final result:
  for (d=0; d<9; ++d) inv[d] = b[d];
}

inline void careful_inverse (const double a[9], double inv[9]) 
{
  double err;
  careful_inverse(a, inv, err);
}
#endif
