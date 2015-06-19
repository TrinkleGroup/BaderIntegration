#ifndef __3DFFT_H
#define __3DFFT_H
#include <math.h>
#include <gsl/gsl_fft_complex.h>

/*
  Program: 3dfft.H
  Author:  D. Trinkle
  Date:    13 July 2005
  Purpose: 3D FFT and iFFT algorithms, where we call GSL for the heavy lifting.
*/

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

// A grid N0 x N1 x N2.  Our indexing of Z is (fortran?) compatible with
// CHGCAR: Z[i + N0*(j + N1*k)] == Z[i][j][k]
// Assume EVERYTHING is already allocated.
int forward_3dfft (double* Z, int N0, int N1, int N2, double* Xr, double* Xi) 
{
  // Wavetables and workspace:
  gsl_fft_complex_wavetable* wave[3];
  wave[0] = gsl_fft_complex_wavetable_alloc(N0);
  wave[1] = gsl_fft_complex_wavetable_alloc(N1);
  wave[2] = gsl_fft_complex_wavetable_alloc(N2);
  
  gsl_fft_complex_workspace* work[3];
  work[0] = gsl_fft_complex_workspace_alloc(N0);
  work[1] = gsl_fft_complex_workspace_alloc(N1);
  work[2] = gsl_fft_complex_workspace_alloc(N2);
  
  // dataspace: ("packed" gsl complex arrays)
  double* data[3], *d;
  data[0] = new double[2*N0];
  data[1] = new double[2*N1];
  data[2] = new double[2*N2];
  
  int i, j, k;
  int nind, N01=N0*N1;
  
  for (k=0; k<N2; ++k) {

    d = data[0];
    for (j=0; j<N1; ++j) {
      // transform Z(i)[j,k] -> Z~(i)[j,k]
      nind = N0*(j+N1*k);
      for (i=0; i<N0; ++i) {
	REAL(d, i) = Z[i+nind];
	IMAG(d, i) = 0;
      }
      gsl_fft_complex_forward(d, 1, N0, wave[0], work[0]);
      for (i=0; i<N0; ++i) {
	Xr[i+nind] = REAL(d, i);
	Xi[i+nind] = IMAG(d, i);
      }
    }

    d = data[1];
    for (i=0; i<N0; ++i) {
      // transform X[i](j)[k] -> X~[i](j)[k]
      nind = i+N0*N1*k;
      for (j=0; j<N1; ++j) {
	REAL(d, j) = Xr[j*N0+nind];
	IMAG(d, j) = Xi[j*N0+nind];
      }
      gsl_fft_complex_forward(d, 1, N1, wave[1], work[1]);
      for (j=0; j<N1; ++j) {
	Xr[j*N0+nind] = REAL(d, j);
	Xi[j*N0+nind] = IMAG(d, j);
      }
    }
  }

  d = data[2];
  for (i=0; i<N0; ++i)
    for (j=0; j<N1; ++j) {
      nind = i + j*N0;
      for (k=0; k<N2; ++k) {
	REAL(d, k) = Xr[k*N01+nind];
	IMAG(d, k) = Xi[k*N01+nind];
      }
      gsl_fft_complex_forward(d, 1, N2, wave[2], work[2]);
      for (k=0; k<N2; ++k) {
	Xr[k*N01+nind] = REAL(d, k);
	Xi[k*N01+nind] = IMAG(d, k);
      }
    }
  
  // Garbage collection:
  for (i=0; i<3; ++i) delete[] data[i];
  for (i=0; i<3; ++i) gsl_fft_complex_workspace_free(work[i]);
  for (i=0; i<3; ++i) gsl_fft_complex_wavetable_free(wave[i]);

  return 0;
}


// Note: this is *destructive* with Xr and Xi.
int inverse_3dfft (double* Xr, double* Xi, int N0, int N1, int N2, double* Z) 
{
  // Wavetables and workspace:
  gsl_fft_complex_wavetable* wave[3];
  wave[0] = gsl_fft_complex_wavetable_alloc(N0);
  wave[1] = gsl_fft_complex_wavetable_alloc(N1);
  wave[2] = gsl_fft_complex_wavetable_alloc(N2);
  
  gsl_fft_complex_workspace* work[3];
  work[0] = gsl_fft_complex_workspace_alloc(N0);
  work[1] = gsl_fft_complex_workspace_alloc(N1);
  work[2] = gsl_fft_complex_workspace_alloc(N2);
  
  // dataspace: ("packed" gsl complex arrays)
  double* data[3], *d;
  data[0] = new double[2*N0];
  data[1] = new double[2*N1];
  data[2] = new double[2*N2];
  
  int i, j, k;
  int nind, N01=N0*N1;
  
  for (k=0; k<N2; ++k) {

    d = data[0];
    for (j=0; j<N1; ++j) {
      // transform Z(i)[j,k] -> Z~(i)[j,k]
      nind = N0*(j+N1*k);
      for (i=0; i<N0; ++i) {
	REAL(d, i) = Xr[i+nind];
	IMAG(d, i) = Xi[i+nind];;
      }
      gsl_fft_complex_inverse(d, 1, N0, wave[0], work[0]);
      for (i=0; i<N0; ++i) {
	Xr[i+nind] = REAL(d, i);
	Xi[i+nind] = IMAG(d, i);
      }
    }

    d = data[1];
    for (i=0; i<N0; ++i) {
      // transform X[i](j)[k] -> X~[i](j)[k]
      nind = i+N0*N1*k;
      for (j=0; j<N1; ++j) {
	REAL(d, j) = Xr[j*N0+nind];
	IMAG(d, j) = Xi[j*N0+nind];
      }
      gsl_fft_complex_inverse(d, 1, N1, wave[1], work[1]);
      for (j=0; j<N1; ++j) {
	Xr[j*N0+nind] = REAL(d, j);
	Xi[j*N0+nind] = IMAG(d, j);
      }
    }
  }

  // Last step--finally overwrite Z:
  int is_real=1;
  double sum_imag = 0;
  d = data[2];
  for (i=0; i<N0; ++i)
    for (j=0; j<N1; ++j) {
      nind = i + j*N0;
      for (k=0; k<N2; ++k) {
	REAL(d, k) = Xr[k*N01+nind];
	IMAG(d, k) = Xi[k*N01+nind];
      }
      gsl_fft_complex_inverse(d, 1, N2, wave[2], work[2]);
      for (k=0; k<N2; ++k) {
	Z[k*N01+nind] = REAL(d, k);
	sum_imag += IMAG(d,k) * IMAG(d,k);
	is_real = is_real && (fabs(IMAG(d, k) < 1e-8));
      }
    }

  sum_imag = sqrt(sum_imag/(N01*N2));
  if (! is_real)
    fprintf(stderr, "## iFFT had rms imaginary component: %.8le\n", sum_imag);
  
  // Garbage collection:
  for (i=0; i<3; ++i) delete[] data[i];
  for (i=0; i<3; ++i) gsl_fft_complex_workspace_free(work[i]);
  for (i=0; i<3; ++i) gsl_fft_complex_wavetable_free(wave[i]);

  return !is_real;
}


#undef REAL
#undef IMAG

#endif
