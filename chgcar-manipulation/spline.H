#ifndef __SPLINE_H
#define __SPLINE_H

/*
  Program: spline.H
  Author:  D. Trinkle
  Date:    August 19, 2004
  Purpose: Define a quartic spline fit to a set of forces along a line.
*/

// ************************** COMPILIATION OPTIONS ************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dcomp.H"
#include "io.H"

// ============================ DATA STRUCTURES ========================

typedef struct 
{
  double a,b,c,d;
} cubic_type;

typedef struct 
{
  double e,a,b,c,d;
} quartic_type;


// ======================== EVALUATION SUBROUTINES =====================
// eval_spline (N, dx, cu_i[], x)  -- evaluates value with x, assuming
//                                    PBC [0, N*dx]

// Take in a u value between 0 and 1, and evaluate:
inline double eval_spline (const cubic_type &cu, const double &u) 
{ return cu.a + u*(cu.b + u*(cu.c + u*cu.d)); }

inline double eval_spline (const quartic_type &qu, const double &u) 
{ return qu.e + u*(qu.a + u*(qu.b + u*(qu.c + u*qu.d))); }


inline void x2ui (const double& x, const double& dx, const int& N, 
		  double& u, int& i) 
{
  u = x/dx;
  i = (int)u;
  u -= i;
  if (u<0) {
    u += 1;
    i -= 1;
  }
  // Now, u should be in [0,1].  We need to put i in [0, N-1]:
  i = i%N;
  if (i<0) i += N;
}

inline double eval_spline (const int &N, const double &dx, cubic_type *cu_i,
			   double x) 
{
  double u;
  int i;
  x2ui(x, dx, N, u, i);
  // Now, pass off to our evaluator:
  return eval_spline(cu_i[i], u);
}

inline double eval_spline (const int &N, const double &dx, quartic_type *qu_i,
			   double x) 
{
  double u;
  int i;
  x2ui(x, dx, N, u, i);
  // Now, pass off to our evaluator:
  return eval_spline(qu_i[i], u);
}


// ========================= MODIFICATION SUBROUTINES ==================
inline void scale_spline(const int &N, cubic_type *cu_i, double scale) 
{
  for (int i=0; i<N; ++i) {
    cu_i[i].a *= scale;
    cu_i[i].b *= scale;
    cu_i[i].c *= scale;
    cu_i[i].d *= scale;
  }
}

inline void scale_spline(const int &N, quartic_type *qu_i, double scale) 
{
  for (int i=0; i<N; ++i) {
    qu_i[i].e *= scale;
    qu_i[i].a *= scale;
    qu_i[i].b *= scale;
    qu_i[i].c *= scale;
    qu_i[i].d *= scale;
  }
}



// ========================= CONSTRUCTION SUBROUTINES ==================

// Constructs the cubic spline for forces, and quartic spline for energies.
// Assumes periodic boundary conditions, and equally spaced points.
// N: number of evaluation points
// dx: spacing between points
// f_i: force at each point
// e0: energy value at 0
//
// e = -df/dx = -df/du * (1/dx)  <- note the sign.
void construct_spline (const int &N, const double &dx, double *f_i, double e0,
		       cubic_type* &f_cu, quartic_type* &e_qu) 
{
  int a, b, k;
  // 1.a. construct the cubic spline propogator to get derivatives of f_i:
  double g[N];
  
  for (a=0; a<N; ++a) g[a] = 0;
  double theta = 2.*M_PI/N;
  for (k=0; k<N; ++k) {
    double pre = -(3./N)*sin(theta*k)/(2+cos(theta*k));
    for (a=0; a<N; ++a)
      g[a] += pre*sin(theta*k*a);
  }
  // 1.b. evaluate derivatives of f_i:
  double fp_i[N];
  for (a=0; a<N; ++a) {
    fp_i[a] = 0;
    for (b=0; b<N; ++b) fp_i[a] += f_i[b] * g[ (a-b+N)%N ];
  }
  
  // 2. fill in values of cubic spline for f, and part of quartic for e:
  f_cu = new cubic_type [N];
  e_qu = new quartic_type [N];
  for (a=0; a<N; ++a) {
    cubic_type* f_p = f_cu + a;
    quartic_type* e_p = e_qu + a;
    f_p->a = f_i[a];
    f_p->b = fp_i[a];
    f_p->c = -3*f_i[a] - 2*fp_i[a] + 3*f_i[(a+1)%N] - fp_i[(a+1)%N];
    f_p->d =  2*f_i[a] + 1*fp_i[a] - 2*f_i[(a+1)%N] + fp_i[(a+1)%N];

    e_p->a =  -dx * f_p->a;
    e_p->b =  -dx * f_p->b / 2.;
    e_p->c =  -dx * f_p->c / 3.;
    e_p->d =  -dx * f_p->d / 4.;
  }
  
  // 3. Propogate forward for e_qu[i].e values.
  e_qu[0].e = e0;
  for (a=0; a<(N-1); ++a)
    e_qu[a+1].e = eval_spline(e_qu[a], 1.);
  
  // 4. Final sanity check: did we stay periodic with e?
  if (! zero( (eval_spline(e_qu[N-1], 1.) - e_qu[0].e) / dx ) )
    fprintf(stderr, "Encountered bad bad problem in construct_spline; should never happen.\n");
}


// just a straight-up energy spline:
void construct_spline (const int &N, const double &dx, double *e_i,
		       cubic_type* &e_cu) 
{
  int a, b, k;
  // 1.a. construct the cubic spline propogator to get derivatives of f_i:
  double g[N];
  
  for (a=0; a<N; ++a) g[a] = 0;
  double theta = 2.*M_PI/N;
  for (k=0; k<N; ++k) {
    double pre = -(3./N)*sin(theta*k)/(2+cos(theta*k));
    for (a=0; a<N; ++a)
      g[a] += pre*sin(theta*k*a);
  }
  // 1.b. evaluate derivatives of e_i:
  double f_i[N];
  for (a=0; a<N; ++a) {
    f_i[a] = 0;
    for (b=0; b<N; ++b) f_i[a] += e_i[b] * g[ (a-b+N)%N ];
  }
  
  // 2. fill in values of cubic spline for e
  if (e_cu == NULL) e_cu = new cubic_type [N];
  for (a=0; a<N; ++a) {
    cubic_type* e_p = e_cu + a;
    e_p->a = e_i[a];
    e_p->b = f_i[a];
    e_p->c = -3*e_i[a] - 2*f_i[a] + 3*e_i[(a+1)%N] - f_i[(a+1)%N];
    e_p->d =  2*e_i[a] + 1*f_i[a] - 2*e_i[(a+1)%N] + f_i[(a+1)%N];
  }
}


// ============================= I/O SUBROUTINES =======================
void write_spline (FILE* file, const int &N, const double &dx, 
		   cubic_type* f_cu, quartic_type* e_qu) 
{
  fprintf(file, "%d %.15le  # number of grid points, spacing\n", N, dx);
  for (int i=0; i<N; ++i) {
    quartic_type* te = e_qu + i;
    fprintf(file, "%23.15le %23.15le %23.15le %23.15le %23.15le\n",
	    te->e, te->a, te->b, te->c, te->d);
    cubic_type* tf = f_cu + i;
    fprintf(file, "%23.15le %23.15le %23.15le %23.15le\n",
	    tf->a, tf->b, tf->c, tf->d);
  }
}

void read_spline (FILE* file, int &N, double &dx, 
		  cubic_type* &f_cu, quartic_type* &e_qu) 
{
  char dump[512];
  nextnoncomment(file, dump, sizeof(dump));
  sscanf(dump, "%d %lf", &N, &dx);
  f_cu = new cubic_type[N];
  e_qu = new quartic_type[N];
  for (int i=0; i<N; ++i) {
    quartic_type* te = e_qu + i;
    nextnoncomment(file, dump, sizeof(dump));
    sscanf(dump, "%lf %lf %lf %lf %lf",
	   &(te->e), &(te->a), &(te->b), &(te->c), &(te->d));
    cubic_type* tf = f_cu + i;
    nextnoncomment(file, dump, sizeof(dump));
    sscanf(dump, "%lf %lf %lf %lf",
	   &(tf->a), &(tf->b), &(tf->c), &(tf->d));
  }
}


#endif
