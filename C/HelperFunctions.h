#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>
#include <fftw3.h>
#include <pthread.h>


/***
*
* --------------------------------------------------------------------------
*
***/
// Dynamically allocate large arrays. C cannot handle large static arrays!!
// double *** = [double**, double**, ...];
// double ** = [double*, double*, ...];
// double * = [double, double, ...];
double **** AllocLargeArray_4D(int Nx, int Ny, int Nz, int depth);
double *** AllocLargeArray_3D(int Nx, int Ny, int Nz);
double ** AllocLargeArray_2D(int Nx, int Ny);

double complex **** AllocLargeArray_4D_Complex(int Nx, int Ny, int Nz, int depth);
double complex *** AllocLargeArray_3D_Complex(int Nx, int Ny, int Nz);
double complex ** AllocLargeArray_2D_Complex(int Nx, int Ny);
/***
*
* --------------------------------------------------------------------------
*
***/



/***
*
* --------------------------------------------------------------------------
*
***/
// Free the dynamically allocated arrays
void Free_4D_mat(int Nx, int Ny, int Nz, int depth, double **** mat);
void Free_3D_mat(int Nx, int Ny, int Nz, double *** mat);
void Free_2D_mat(int Nx, int Ny, double ** mat);

void Free_4D_mat_Complex(int Nx, int Ny, int Nz, int depth, double complex **** mat);
void Free_3D_mat_Complex(int Nx, int Ny, int Nz, double complex *** mat);
void Free_2D_mat_Complex(int Nx, int Ny, double complex ** mat);
/***
*
* --------------------------------------------------------------------------
*
***/



/***
*
* --------------------------------------------------------------------------
*
***/
// Convert 3D array into row major 1D array
void Convert3D_2_RowMajor(int Nx, int Ny, int Nz, double *** in, double * out);
void Convert2D_2_RowMajor(int Nx, int Ny, double ** in, double * out);

// Convert 2D array into col major 1D array
void Convert2D_2_ColMajor(int Nx, int Ny, double ** in, double * out);
void Convert2D_2_ColMajor_Complex(int Nx, int Ny, double complex ** in, double complex * out);
void Convert2D_2_ColMajor_Double_2_Complex(int Nx, int Ny, 
	double ** in, double complex * out);
void Convert3D_2_RowMajor_Complex(int Nx, int Ny, int Nz, double complex *** in, double complex * out);
/***
*
* --------------------------------------------------------------------------
*
***/



/***
*
* --------------------------------------------------------------------------
*
***/
// Convert row major array into 3D
void ConvertRowMajor_2_3D(int Nx, int Ny, int Nz, double * in, double *** out);
/***
*
* --------------------------------------------------------------------------
*
***/


/***
*
* --------------------------------------------------------------------------
*
***/
double CalcError(double complex * prev, double complex * next, int total);

// Print a 3D Matrix out, MATLAB style printing
void Print3DMat(int Nx, int Ny, int Nz, double *** in);
void SaveFile(char * filename, int total, double * array_RowMaj);

// flag = 1, save real and imag
// flag = 0, save real only
void SaveFile_Complex(char * filename, int total, int flag, double complex * array_RowMaj);
/***
*
* --------------------------------------------------------------------------
*
***/

/***
*
* --------------------------------------------------------------------------
*
***/
// meshgrid
void Meshgrid(int Nx, int Ny, int Nz, double * x, double * y, double * z, 
				double *** x_grid, double *** y_grid, double *** z_grid);


double MaxValue(double myArray[], size_t size);
double complex MaxValue_Complex(double complex myArray[], size_t size);
/***
*
* --------------------------------------------------------------------------
*
***/



/***
*
* --------------------------------------------------------------------------
*
***/
// LAPACK wrappers
// Require ColMaj matrices
void InverseMat(double * A, int N);
void InverseMat_Complex(double complex * A, int N);
void MatrixComplexEigensystem( double complex * eigenvectorsVR, double complex * eigenvaluesW, double complex * A, int N);
void LU_Decomp_Complex(int size, double complex * A, double complex * L,
				double complex * U, double * P);
/***
*
* --------------------------------------------------------------------------
*
***/



/***
*
* --------------------------------------------------------------------------
*
***/
void DotMult(int N, double * A, double * B, double * C);
void DotMult_Complex(int N, double complex * A, double complex * B, double complex * C);
/***
*
* --------------------------------------------------------------------------
*
***/



/***
*
* --------------------------------------------------------------------------
*
***/
// FFT wrappers, require Row Major matrices
void FFT3(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array);
void IFFT3(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array);

void FFT2(int Nx, int Ny, 
	double complex * in_array, double complex * out_array);
void IFFT2(int Nx, int Ny,
	double complex * in_array, double complex * out_array);
	
void IFFT2_slices(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array, int start, int end);
void IFFT2_slices_MultiThreaded(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array, int NUM_THREADS);
/***
*
* --------------------------------------------------------------------------
*
***/
	
#endif