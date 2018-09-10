#ifndef __LINEAR_H
#define __LINEAR_H
#include <math.h>

/*
  Program: linear.H
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: Linear interpolation algorithm on a grid
*/

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])


const int TRANSLATION[8][3] = {
  {0,0,0},{0,0,1},{0,1,0},{0,1,1},
  {1,0,0},{1,0,1},{1,1,0},{1,1,1} };

const double DIFF_TRANSLATION[8][3] = {
  {1,1,1},{1,1,0},{1,0,1},{1,0,0},
  {0,1,1},{0,1,0},{0,0,1},{0,0,0} };
    

// Z is defined on a grid, given by N:
// Z[n0 + N[0]*(n1 + N[1]*n2)]
double linear (double* Z, int N[3], double r[3]) 
{
  int i, n[3];
  double x[3];
  for (i=0; i<3; ++i) {
    x[i]=r[i];
    while (x[i]<0) x[i] += N[i];
    while (x[i]>=N[i]) x[i] -= N[i];
    n[i] = (int)(x[i]);
    x[i] = x[i] - n[i];
  }
  // Now, n[i] should be >=0 and <N[i] and
  //      x[i] should be >=0 and <1
  
  int index[8][3];
  double scale[8];
  for (int t=0; t<8; ++t) {
    for (i=0; i<3; ++i)
      index[t][i] = (n[i] + TRANSLATION[t][i])%(N[i]);
    scale[t] = fabs( (DIFF_TRANSLATION[t][0]-x[0])*(DIFF_TRANSLATION[t][1]-x[1])*
		     (DIFF_TRANSLATION[t][2]-x[2]) );
  }
  double sum=0;
  int N0=N[0], N01=N[0]*N[1];
  for (int t=0; t<8; ++t)
    if (scale[t]>1e-10)
      sum += scale[t]*Z[index[t][0] + N0*index[t][1] + N01*index[t][2]];
  return sum;
}

#endif
