#include "fftw3/include/fftw3.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Convert 3d array into row major array
void convert3D_2_RowMajor(int Nx, int Ny, int Nz, double *** in, double * out)
{
	for(int i = 0; i < Nx, i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		int index = k + Nz * (j + Ny * i);
		out[index] = in[i][j][k];
	}}}
}

void convertRowMajor_2_3D(int Nx, int Ny, int Nz, double * in, double *** out)
{
	for(int i = 0; i < Nx, i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		int index = k + Nz * (j + Ny * i);
		out[i][j][k] = in[index];
	}}}
}

int main()
{
	// Test 3D
	// row major
	int Nx = 4; int Ny = 4; int Nz = 4;
		
	


	int N = 1e8;
	fftw_complex *in, *out;
	fftw_plan fft;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	for(int i = 0; i < N; i++)
	{
		in[i][0] = i;
		in[i][1] = i * 3;
//		fprintf(stdout, "%g %gi\n", in[i][0], in[i][1]);
	}	

//	fprintf(stdout, "=======================\n");

	fft = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(fft); /* repeat as needed */

	for(int i = 0; i < N; i++)
	{
//		fprintf(stdout, "%g %gi\n", out[i][0], out[i][1]);
	}


//	fprintf(stdout, "=======================\n");

	double scale = N;

	fftw_plan ifft;
	ifft = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(ifft);

	for(int i = 0; i < N; i++)
	{
//		fprintf(stdout, "%g %gi\n", in[i][0]/scale, in[i][1]/scale);
	}

	fftw_destroy_plan(fft);
	fftw_free(in); fftw_free(out);
}
