#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fftw3/include/fftw3.h"
#include <string.h>

#include "Setup.h"

// Dynamically allocate large arrays. C cannot handle large static arrays!!
// double *** = [double**, double**, ...];
// double ** = [double*, double*, ...];
// double * = [double, double, ...];
double *** AllocLargeArray(double Nx, double Ny, double Nz)
{
	double *** out = (double ***) malloc(sizeof(double **) * Nx);
	for(int i = 0; i < Nx; i ++)
	{
		out[i] = (double **) malloc(sizeof(double *) * Ny);
		for(int j = 0; j < Ny; j++)
		{
			out[i][j] = (double *) malloc(sizeof(double) * Nz);
			for(int k = 0; k < Nz; k++)
			{
				out[i][j][k] = 0; // Initialize everything to zero
			}
		}
	}
	return out;
}

// Convert 3d array into row major array
void convert3D_2_RowMajor(int Nx, int Ny, int Nz, double *** in, double * out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		int index = k + Nz * (j + Ny * i);
		out[index] = in[i][j][k];
	}}}
}

// Convert row major array into 3D
void convertRowMajor_2_3D(int Nx, int Ny, int Nz, double * in, double *** out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		int index = k + Nz * (j + Ny * i);
		out[i][j][k] = in[index];
	}}}
}


void print3DMat(int Nx, int Ny, int Nz, double *** in)
{
	for(int k = 0; k < Nz; k++){
		fprintf(stdout, "Depth %u\n", k);
		for(int i = 0; i < Nx; i++){
			for(int j = 0; j < Ny; j++){
				fprintf(stdout, "%g ", in[i][j][k]);
			}
			fprintf(stdout, "\n");
		}
	}
}


// Main function
int main()
{
	// Set up randomization
	srand(1);
	int r = rand();
	
	// Set up Constants
	struct ConstantsStruct Constants;
	InitSetup(&Constants);

	unsigned int Nx = Constants.Nx;
	unsigned int Ny = Constants.Ny;
	unsigned int Nz = Constants.Nz;
	
	double * x_axis = Constants.x_axis;
	double * y_axis = Constants.y_axis;
	double * z_axis = Constants.z_axis;
	
	// --- Set up Axes --- //
	double *** x_grid = AllocLargeArray(Nx, Ny, Nz);
 	double *** y_grid = AllocLargeArray(Nx, Ny, Nz);
	double *** z_grid = AllocLargeArray(Nx, Ny, Nz);

 	// Create x, y, and z grids
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		x_grid[i][j][k] = x_axis[i];
	}}}

 	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		y_grid[i][j][k] = y_axis[j];
	}}}

	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		z_grid[i][j][k] = z_axis[k];
	}}}
	
	// 
	print3DMat(Nx, Ny, Nz, x_grid);
	
	// 3D FFT Test! 
	unsigned int Total = Nx * Ny * Nz;

	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Total);
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Total);
	double * in_array = malloc(sizeof(double) * Total);
	convert3D_2_RowMajor(Nx, Ny, Nz, x_grid, in_array);
	for(int i = 0; i < Total; i++)
	{
		in[i][0] = in_array[i];
		in[i][1] = 0;
	}
	fftw_plan FFT_3D = fftw_plan_dft_3d(Nx, Ny, Nz,
                           in, out,
                           FFTW_FORWARD, FFTW_ESTIMATE);
						   
	fftw_execute(FFT_3D);
	
	double * out_array_real = malloc(sizeof(double) * Total);
	double * out_array_imag = malloc(sizeof(double) * Total);

	for(int i = 0; i < Total; i++)
	{
		out_array_real[i] = out[i][0];
		out_array_imag[i] = out[i][1];
	}
	
	double *** out_real_3D = AllocLargeArray(Nx, Ny, Nz);
	double *** out_imag_3D = AllocLargeArray(Nx, Ny, Nz);
	
	convertRowMajor_2_3D(Nx, Ny, Nz, out_array_real, out_real_3D);
	convertRowMajor_2_3D(Nx, Ny, Nz, out_array_imag, out_imag_3D);
	
	
	fprintf(stdout, "============\n");
	print3DMat(Nx, Ny, Nz, out_real_3D);
	fprintf(stdout, "============\n");
	print3DMat(Nx, Ny, Nz, out_imag_3D);
	
	// 3D IFFT Test
	fftw_plan IFFT_3D = fftw_plan_dft_3d(Nx, Ny, Nz,
                           out, in,
                           FFTW_BACKWARD, FFTW_ESTIMATE);
						   
	fftw_execute(IFFT_3D);	
	for(int i = 0; i < Total; i++)
	{
		out_array_real[i] = in[i][0] / Total;
		out_array_imag[i] = in[i][1] / Total;
	}
	
	convertRowMajor_2_3D(Nx, Ny, Nz, out_array_real, out_real_3D);
	convertRowMajor_2_3D(Nx, Ny, Nz, out_array_imag, out_imag_3D);
	
	fprintf(stdout, "============\n");
	print3DMat(Nx, Ny, Nz, out_real_3D);
	fprintf(stdout, "============\n");
	print3DMat(Nx, Ny, Nz, out_imag_3D);
		
	
	double complex * test = calloc(3 * 3 * 3, sizeof(double complex));
	double complex * out = calloc(27, sizeof(double complex));
	for(int i = 0; i < 27; i++)
		test[i] = i;
	
	FFT3(3, 3, 3, test, out);
	
	double * out_array_real = malloc(sizeof(double) * 27);

	for(int i = 0; i < 27; i++)
	{
		printf("%g + %gi\n", creal(out[i]), cimag(out[i]));
		out_array_real[i] = creal(out[i]);
	}
	double *** out_real_3D = AllocLargeArray_3D(3,3,3);
	
	ConvertRowMajor_2_3D(3, 3, 3, out_array_real, out_real_3D);
	Print3DMat(3, 3, 3, out_real_3D);	
	return 0;
}