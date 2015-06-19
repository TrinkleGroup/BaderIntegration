/*
  Header:  ws_voronoi.H
  Author:  D. Trinkle
  Date:    2010 December 27
  Purpose: Determines the prefactors for computation of flux in Wigner-Seitz
           grid cells, based on the Voronoi decomposition.
*/

#include <stdlib.h>

const double TOLER=1e-8;

// is a vector r inside of the WS cell defined by the set of vectors R?
// NOTE: R[n][0..2] = cartesian coord, R[n][3] = (R.R)/2
inline int incell (double r[3], int Nneigh, double** R) 
{
  for (int n=0; n<Nneigh; ++n) 
    if (dot(r, R[n]) > (R[n][3]+TOLER)) return 0;
  return 1;
}

// want this list in ascending order 
int vert_comp(const void* a, const void* b) 
{
  double *av=*((double**)a);
  double *bv=*((double**)b);
  if ( av[3] < bv[3] ) return -1;
  if ( av[3] > bv[3] ) return  1;
  return 0;
}

// Construct a list of the neighboring vectors that define the Wigner-Seitz cell
// and compute the "alpha" factors needed for flux; you multiply the difference
// in densities by alpha, and use this to compute the transition probabilities.
int gen_WS_voronoi(double alatt[9], int &Nvect, int** &vect, double* &alpha) 
{
  // 0. Generate (pruned) list of neighboring vectors that bound the
  // Wigner-Seitz cell.
  // Note for future: should precompute this by making a sphere of radius
  // with the largest length vector multiplied by, say, 2.
  const int Nrange=3;	// (generally) overkill; need a way to compute this
  int Nneigh = (2*Nrange+1)*(2*Nrange+1)*(2*Nrange+1)-1;
  double **R = new double*[Nneigh];
  int **nvect = new int*[Nneigh];
  Nneigh=0;
  { int nv[3];
    for (nv[0]=-Nrange; nv[0]<=Nrange; ++nv[0])
      for (nv[1]=-Nrange; nv[1]<=Nrange; ++nv[1])
	for (nv[2]=-Nrange; nv[2]<=Nrange; ++nv[2])
	  if (! zero_vect(nv) ) {
	    nvect[Nneigh]=new int[3];
	    R[Nneigh]=new double[4];
	    for (int d=0; d<3; ++d) nvect[Nneigh][d]=nv[d];
	    mult_vect(alatt, nv, R[Nneigh]);
	    R[Nneigh][3] = 0.5*dot(R[Nneigh], R[Nneigh]);
	    Nneigh++;
	  }
  }
  // prune that list: basically, if R/2 isn't inside the WS cell, then
  // R isn't involved in bounding the WS cell.
  for (int n=(Nneigh-1); n>=0; --n) {
    double r[3] = {0.5*R[n][0], 0.5*R[n][1], 0.5*R[n][2]};
    if (! incell(r, Nneigh, R) ) {
      // need to remove this entry...
      delete[] nvect[n];
      delete[] R[n];
      for (int i=n; i<(Nneigh-1); ++i) {
	nvect[i] = nvect[i+1];
	R[i] = R[i+1];
      }
      --Nneigh;
    }
  }
  
  // 1. Run over each facet to find all of the vertex points
  int MAXVERT = (Nneigh-2)*(Nneigh-4);
  double** rvert = new double*[MAXVERT];
  for (int nv=0; nv<MAXVERT; ++nv) rvert[nv] = new double[4]; // [3] = phi
  double* alph = new double[Nneigh];
  Nvect=Nneigh;
  for (int n=0; n<Nneigh; ++n) {
    int nvert = 0;
    double Rdot[9], detR, Rinv[9], R2[3];
    for (int d=0; d<3; ++d) Rdot[d] = R[n][d];
    R2[0] = R[n][3];
    for (int nA=0; nA<Nneigh; ++nA) {
      for (int d=0; d<3; ++d) Rdot[3+d] = R[nA][d];
      R2[1] = R[nA][3];
      for (int nB=(nA+1); nB<Nneigh; ++nB) {
	for (int d=0; d<3; ++d) Rdot[6+d] = R[nB][d];
	R2[2] = R[nB][3];
	detR = inverse(Rdot, Rinv);
	if (fabs(detR)>TOLER) {
	  mult_vect(Rinv, R2, rvert[nvert]);
	  for (int d=0; d<3; ++d) rvert[nvert][d] *= 1./detR;
	  if (incell(rvert[nvert], Nneigh, R) ) ++nvert;
	}
      }
    }
    // check to make sure none of the vertices correspond to R/2:
    int zeroarea=0;
    for (int nv=0; nv<nvert; ++nv)
      if (fabs(dot(rvert[nv], rvert[nv])-0.5*R[n][3])<TOLER) zeroarea=1;
    if (zeroarea || (nvert==0)) {
      alph[n]=0;
      --Nvect;
      continue;
    }
    // Now we have a list of all the vertices for the polygon
    // defining the facet along the direction R[n].
    // Last step is to sort the list in terms of a winding angle around
    // R[n].  To do that, we define rx and ry which are perpendicular
    // to R[n], normalized, and right-handed: ry = R x rx, so that
    // rx x ry points along R[n].
    double rx[3], ry[3];
    for (int d=0; d<3; ++d) rx[d] = rvert[0][d];
    double rdRn = dot(rx, R[n])/dot(R[n],R[n]);
    for (int d=0; d<3; ++d) rx[d] -= rdRn*R[n][d];
    rdRn = sqrt(dot(rx, rx));
    for (int d=0; d<3; ++d) rx[d] *= 1./rdRn;
    crossprod(R[n], rx, ry);
    rdRn = sqrt(dot(ry, ry));
    for (int d=0; d<3; ++d) ry[d] *= 1./rdRn;
    // now compute winding angle phi_n
    for (int nv=0; nv<nvert; ++nv)
      rvert[nv][3] = atan2(dot(rvert[nv], ry), dot(rvert[nv], rx));
    // sort that list
    qsort(rvert, nvert, sizeof(double*), vert_comp);
    alph[n] = 0;
    for (int nv=0; nv<nvert; ++nv)
      alph[n] += tripleprod(rvert[nv], rvert[(nv+1)%nvert], R[n]);
    alph[n] *= 0.25/R[n][3];
    if (fabs(alph[n])<TOLER) {
      alph[n]=0;
      --Nvect;
    }
  }

  // 2. Output
  vect = new int*[Nvect];
  alpha = new double[Nvect];
  int nv=0;
  for (int n=0; n<Nneigh; ++n)
    if (alph[n]!=0) {
      vect[nv] = new int[3];
      for (int d=0; d<3; ++d) vect[nv][d] = nvect[n][d];
      alpha[nv] = alph[n];
      ++nv;
    }
  
  // 3. Garbage collection
  for (int nv=0; nv<MAXVERT; ++nv) delete[] rvert[nv];
  delete[] rvert;
  for (int n=0; n<Nneigh; ++n) {
    delete[] nvect[n];
    delete[] R[n];
  }
  delete[] nvect;
  delete[] R;
  delete[] alph;
  
  return 0;
}
