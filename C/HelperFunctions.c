#include "HelperFunctions.h"

double CalcError(double complex * prev, double complex * next, int total)
{
	double out = 0;
	for(int i = 0; i < total; i++)
	{
		out = out + fabs(creal(prev[i]) - creal(next[i]));
	}
	
	return out / (double)total;
}

void Print1DVec(double * in, int size)
{
	for(int i = 0; i < size; i++)
		fprintf(stdout, "%.4g\t", in[i]);
	fprintf(stdout, "\n");
}

// flag = 1, save real and imag
// flag = 0, save real only
void SaveFile_Complex(char * filename, int total, int flag, double complex * array_RowMaj)
{
	if(flag == 1)
	{
		FILE * out = fopen(filename, "w");
		for(int i = 0; i < total; i++)
			fprintf(out, "%g + %g * 1i\n", creal(array_RowMaj[i]), cimag(array_RowMaj[i]));	
		fclose(out);
	}
	else
	{
		FILE * out = fopen(filename, "w");
		for(int i = 0; i < total; i++)
			fprintf(out, "%g\n", creal(array_RowMaj[i]));	
		fclose(out);	
	}	
}

void SaveFile(char * filename, int total, double * array_RowMaj)
{
	FILE * out = fopen(filename, "w");
	for(int i = 0; i < total; i++)
		fprintf(out, "%g\n", array_RowMaj[i]);	
	fclose(out);
}

/***
*
* --------------------------------------------------------------------------
*
***/
// Dynamically allocate large arrays. C cannot handle large static arrays!!
// double *** = [double**, double**, ...];
// double ** = [double*, double*, ...];
// double * = [double, double, ...];
double **** AllocLargeArray_4D(int Nx, int Ny, int Nz, int depth)
{
	double **** out = malloc(sizeof(double ***) * Nx);
	for(int i = 0; i < Nx; i ++)
	{
		out[i] = malloc(sizeof(double **) * Ny);
		for(int j = 0; j < Ny; j++)
		{
			out[i][j] = malloc(sizeof(double *) * Nz);
			for(int k = 0; k < Nz; k++)
			{
				out[i][j][k] = malloc(sizeof(double) * depth);
				for(int l = 0; l < depth; l++)
				{
					out[i][j][k][l] = 0; // Initialize everything to zero
				}
			}
		}
	}
	return out;
}

double *** AllocLargeArray_3D(int Nx, int Ny, int Nz)
{
	double *** out = malloc(sizeof(double **) * Nx);
	for(int i = 0; i < Nx; i ++)
	{
		out[i] = malloc(sizeof(double *) * Ny);
		for(int j = 0; j < Ny; j++)
		{
			out[i][j] = malloc(sizeof(double) * Nz);
			for(int k = 0; k < Nz; k++)
			{
				out[i][j][k] = 0; // Initialize everything to zero
			}
		}
	}
	return out;
}

double ** AllocLargeArray_2D(int Nx, int Ny)
{
	double ** out = malloc(sizeof(double *) * Nx);
	for(int i = 0; i < Nx; i++)
	{
		out[i] = malloc(sizeof(double) * Ny);
		for(int j = 0; j < Ny; j++)
		{
			out[i][j] = 0;
		}
	}
	return out;
}

double complex **** AllocLargeArray_4D_Complex(int Nx, int Ny, int Nz, int depth)
{
	double complex **** out = malloc(sizeof(double complex ***) * Nx);
	for(int i = 0; i < Nx; i ++)
	{
		out[i] = malloc(sizeof(double complex **) * Ny);
		for(int j = 0; j < Ny; j++)
		{
			out[i][j] = malloc(sizeof(double complex *) * Nz);
			for(int k = 0; k < Nz; k++)
			{
				out[i][j][k] = malloc(sizeof(double complex) * depth);
				for(int l = 0; l < depth; l++)
				{
					out[i][j][k][l] = 0; // Initialize everything to zero
				}
			}
		}
	}
	return out;
}

double complex *** AllocLargeArray_3D_Complex(int Nx, int Ny, int Nz)
{
	double complex *** out = malloc(sizeof(double complex **) * Nx);
	for(int i = 0; i < Nx; i ++)
	{
		out[i] = malloc(sizeof(double complex *) * Ny);
		for(int j = 0; j < Ny; j++)
		{
			out[i][j] = malloc(sizeof(double complex) * Nz);
			for(int k = 0; k < Nz; k++)
			{
				out[i][j][k] = 0; // Initialize everything to zero
			}
		}
	}
	return out;
}

double complex ** AllocLargeArray_2D_Complex(int Nx, int Ny)
{
	double complex ** out = malloc(sizeof(double complex *) * Nx);
	for(int i = 0; i < Nx; i++)
	{
		out[i] = malloc(sizeof(double complex) * Ny);
		for(int j = 0; j < Ny; j++)
		{
			out[i][j] = 0;
		}
	}
	return out;
}

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
void Free_4D_mat(int Nx, int Ny, int Nz, int depth, double **** mat)
{
	for(int i = 0; i < Nx; i ++)
	{
		Free_3D_mat(Ny, Nz, depth, mat[i]);
	}
	free(mat);
}

void Free_3D_mat(int Nx, int Ny, int Nz, double *** mat)
{
	for(int i = 0; i < Nx; i++)
	{
		Free_2D_mat(Ny, Nz, mat[i]);
	}
	free(mat);
}

void Free_2D_mat(int Nx, int Ny, double ** mat)
{
	for(int i = 0; i < Nx; i++){
		free(mat[i]);
	}
	free(mat);
}

void Free_4D_mat_Complex(int Nx, int Ny, int Nz, int depth, double complex **** mat)
{
	for(int i = 0; i < Nx; i ++)
	{
		Free_3D_mat_Complex(Ny, Nz, depth, mat[i]);
	}
	free(mat);
}

void Free_3D_mat_Complex(int Nx, int Ny, int Nz, double complex *** mat)
{
	for(int i = 0; i < Nx; i++)
	{
		Free_2D_mat_Complex(Ny, Nz, mat[i]);
	}
	free(mat);
}

void Free_2D_mat_Complex(int Nx, int Ny, double complex ** mat)
{
	for(int i = 0; i < Nx; i++){
		free(mat[i]);
	}
	free(mat);
}
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
// Convert 3d array into row major array
void Convert3D_2_RowMajor(int Nx, int Ny, int Nz, double *** in, double * out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		int index = k + Nz * (j + Ny * i);
		out[index] = in[i][j][k];
	}}}
}

void Convert2D_2_RowMajor(int Nx, int Ny, double ** in, double * out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		out[j + Ny * i] = in[i][j];
	}}
}

void Convert2D_2_ColMajor(int Nx, int Ny, double ** in, double * out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		out[i + Nx * j] = in[i][j];
	}}
}

void Convert2D_2_ColMajor_Complex(int Nx, int Ny, double complex ** in, double complex * out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		out[i + Nx * j] = in[i][j];
	}}
}

void Convert2D_2_ColMajor_Double_2_Complex(int Nx, int Ny, 
	double ** in, double complex * out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		out[i + Nx * j] = in[i][j];
	}}
}
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
void ConvertRowMajor_2_3D(int Nx, int Ny, int Nz, double * in, double *** out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		out[i][j][k] = in[k + Nz * (j + Ny * i)];
	}}}
}

void ConvertRowMajor_2_3D_Complex(int Nx, int Ny, int Nz, double complex * in, double complex *** out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		out[i][j][k] = in[k + Nz * (j + Ny * i)];
	}}}
}

void Convert3D_2_RowMajor_Complex(int Nx, int Ny, int Nz, 
	double complex *** in, double complex * out)
{
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		int index = k + Nz * (j + Ny * i);
		out[index] = in[i][j][k];
	}}}
}
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
// Print a 3D Matrix out, MATLAB style printing
void Print3DMat(int Nx, int Ny, int Nz, double *** in)
{
	for(int k = 0; k < Nz; k++){
		fprintf(stdout, "Depth %u\n", k);
		for(int i = 0; i < Nx; i++){
			for(int j = 0; j < Ny; j++){
				fprintf(stdout, "%.4g\t", in[i][j][k]);
			}
			fprintf(stdout, "\n");
		}
	}
}
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
void Meshgrid(int Nx, int Ny, int Nz, double * x, double * y, double * z, 
				double *** x_grid, double *** y_grid, double *** z_grid)
{
	 // Create x, y, and z grids
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			for(int k = 0; k < Nz; k++){
				x_grid[i][j][k] = x[i];
				y_grid[i][j][k] = y[j];
				z_grid[i][j][k] = z[k];
			}
		}
	}
	
}

double MaxValue(double myArray[], size_t size) 
{
    /* enforce the contract */
    assert(myArray && size);
    size_t i;
    double maxValue = myArray[0];

    for (i = 1; i < size; ++i) {
        if ( myArray[i] > maxValue ) {
            maxValue = myArray[i];
        }
    }
    return maxValue;
}

double Max2(double a, double b)
{
	if( a > b ) return a;
	return b;
}

double complex MaxValue_Complex(double complex myArray[], size_t size) 
{
    /* enforce the contract */
    assert(myArray && size);
    size_t i;
    double complex maxValue = myArray[0];

    for (i = 1; i < size; ++i) {
        if ( Max2(creal(myArray[i]), cimag(myArray[i])) > 
			Max2(creal(maxValue), cimag(maxValue)) )
		{
            maxValue = myArray[i];
        }
    }
    return maxValue;
}

// LU decomposition on square matrices
void LU_Decomp_Complex(int size, double complex * A, double complex * L,
				double complex * U, double * P)
{
	int piv[size];
	int info;
	
	// LU decomp -> column major data
	zgetrf_(&size, &size, A, &size, piv, &info);
	
	if(info < 0)
	{
		printf("Error %d\n", info);
		return;
	}
	
	// Extract L, U
	for(int i = 0; i < size; i++){
	for(int j = 0; j < size; j++){
		if(i > j) // Off diagonal, L
			L[i + j * size] = A[i + j * size];
			
		else // U
			U[i + j * size] = A[i + j * size];
	}}
	for(int i = 0; i < size; i++)
		L[size * i + i] = 1;	
	
	
	
	
	// Extract P
 	
	// start = identity matrix
	double * start = calloc(size * size, sizeof(double));
	for(int i = 0; i < size; i++)
		start[size * i + i] = 1;
	
	double * start_temp = calloc(size * size, sizeof(double));
	memcpy(start_temp, start, size * size * sizeof(double));
	
	double * m = calloc(size * size, sizeof(double));
	double * t1 = calloc(size, sizeof(double));
	double * t2 = calloc(size, sizeof(double));
	
	int n; char N = 'N'; 
	double zero = 0.0;
	double one = 1.0;
		
	// i -> index piv
	// i = 1 ... size
	for(int i = 1; i < size + 1; i++)
	{
		memset(m, 0, sizeof(double) * size * size); // m = identity matrix
		for(int j = 0; j < size; j++)
			m[size * j + j] = 1;
		
		memset(t1, 0, size * sizeof(double));
		memset(t2, 0, size * sizeof(double));
		
		n = piv[i - 1]; // n = 1 ... size
		
		if(n != i)
		{
			// To swap...
			t1[n - 1] = 1;
			t2[i - 1] = 1;
			
			// swap
			for(int j = 0; j < size; j++)
			{
				m[i - 1 + size * j] = t1[j];
				m[n - 1 + size * j] = t2[j];
			}
			
			// start = start * m
			dgemm_(&N, &N, &size, &size, &size, 
					&one, start, &size, 
					m, &size, 
					&zero, start_temp, &size);
			memcpy(start, start_temp, size * size * sizeof(double));
		}
	}
	
	for(int i = 0; i < size; i++){
	for(int j = 0; j < size; j++){
		P[i + j * size] = start[j + i * size];
	}}
	
	free(m); free(start); 
	free(t1); free(t2);
	free(start_temp);

}

/* a = [1 5 3 6 6 6];
start = eye(6);

for i = 1 : 6
    m = eye(6);

    t1 = [];
    t2 = [];
    
    n = a(i);
    if n ~= i
        for j = 1 : n - 1
            t1(j) = 0;
        end
        t1(n) = 1;
        for j = n + 1 : 6
            t1(j) = 0;
        end
        
        for j = 1 : i - 1
            t2(j) = 0;
        end
        t2(i) = 1;
        for j = i + 1 : 6
            t2(j) = 0;
        end
        m(i,:) = t1;
        m(n,:) = t2;
    end
    
    start = start * m;
    
end

start = start'; */

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
//
// LAPACK operates on column major 1D arrays!!!
//
// Find inverse of a 2D array
// A is column major 1D
void InverseMat(double * A, int N)
{
    int *IPIV = malloc(sizeof(int) * (N+1));
    int LWORK = N*N;
    double *WORK = malloc(sizeof(double) * LWORK);
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}

void InverseMat_Complex(double complex * A, int N)
{
    int *IPIV = malloc(sizeof(int) * (N+1));
    int LWORK = N*N;
    double *WORK = malloc(sizeof(double complex) * LWORK);
    int INFO;

    zgetrf_(&N,&N,A,&N,IPIV,&INFO);
    zgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}

// Find eigenvectors and eigenvalues of a complex system
void MatrixComplexEigensystem( double complex * eigenvectorsVR, double complex * eigenvaluesW, 
								double complex * A, int N)
{
	char JOBVL ='N';   // Compute Right eigenvectors
	char JOBVR ='V';   // Do not compute Left eigenvectors
	double complex VL[1];

	int LDVL = 1; 
	int LDVR = N;
	int LWORK = 4*N; 

	double complex *WORK =  (double complex*)malloc( LWORK*sizeof(double complex));
	double complex *RWORK = (double complex*)malloc( 2*N*sizeof(double complex));

	int INFO;

	zgeev_( &JOBVL, &JOBVR, &N, A ,  &N , eigenvaluesW ,

	VL, &LDVL,
	eigenvectorsVR, &LDVR, 
	WORK, 
	&LWORK, RWORK, &INFO );

	free(WORK);
	free(RWORK);
}
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
void DotMult(int N, double * A, double * B, double * C)
{
	for(int i = 0; i < N; i++)
	{
		C[i] = A[i] * B[i];
	}
}
void DotMult_Complex(int N, double complex * A, double complex * B, double complex * C)
{
	for(int i = 0; i < N; i++)
	{
		C[i] = A[i] * B[i];
	}	
}
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
void FFT3(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array)
{
	int total = Nx * Ny * Nz;
	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);

	// Copy in_array into in...
	memcpy(in, in_array, sizeof(fftw_complex) * total);

	fftw_plan FFT_3D = fftw_plan_dft_3d(Nx, Ny, Nz,
                           in, out,
                           FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(FFT_3D);
	
	// Copy into out_array...
	memcpy(out_array, out, sizeof(fftw_complex) * total);
	
	fftw_free(in); fftw_free(out);
	fftw_destroy_plan(FFT_3D);
}

// IFFT
// 3D
// Same as MATLAB, scale by total
// Force to be real
void IFFT3(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array)
{
	int total = Nx * Ny * Nz;
	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);

	// Copy in_array into in...
	memcpy(in, in_array, sizeof(fftw_complex) * total);

	fftw_plan FFT_3D = fftw_plan_dft_3d(Nx, Ny, Nz,
                           in, out,
                           FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(FFT_3D);
	
	// Copy into out_array...
	for(int i = 0; i < total; i++)
	{
		out_array[i] = creal(out[i]) / total;
	}
	
	fftw_free(in); fftw_free(out);
	fftw_destroy_plan(FFT_3D);
}


void FFT2(int Nx, int Ny, 
	double complex * in_array, double complex * out_array)
{
	int total = Nx * Ny;
	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);

	// Copy in_array into in...
	memcpy(in, in_array, sizeof(fftw_complex) * total);

	fftw_plan FFT_2D = fftw_plan_dft_2d(Nx, Ny,
                           in, out,
                           FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(FFT_2D);
	
	// Copy into out_array...
	memcpy(out_array, out, sizeof(fftw_complex) * total);
	
	fftw_free(in); fftw_free(out);
	fftw_destroy_plan(FFT_2D);
}

void IFFT2(int Nx, int Ny,
	double complex * in_array, double complex * out_array)
{
	int total = Nx * Ny;
	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * total);

	// Copy in_array into in...
	memcpy(in, in_array, sizeof(fftw_complex) * total);

	fftw_plan FFT_2D = fftw_plan_dft_2d(Nx, Ny,
                           in, out,
                           FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(FFT_2D);
	
	// Copy into out_array...
	for(int i = 0; i < total; i++)
	{
		out_array[i] = creal(out[i]) / total;
	}
	
	fftw_free(in); fftw_free(out);
	fftw_destroy_plan(FFT_2D);
}


void IFFT2_slices_MultiThread(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array, int start, int end)
{
	
}
	
// IFFT2 each slice of a 3D row major matrix
void IFFT2_slices(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array, int start, int end)
{	
	double complex * temp_in = calloc(Nx * Ny, sizeof(double complex));
	double complex * temp_out = calloc(Nx * Ny, sizeof(double complex));
	for(int z = start; z < end; z++)
	{
		// Transfer in from 3D row major
		for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			temp_in[j + Ny * i] = in_array[z + Nz * (j + Ny * i)];
		}}

		// FFT2
		IFFT2(Nx, Ny, temp_in, temp_out);

		// Transfer out into 3D row major
		for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			out_array[z + Nz * (j + Ny * i)] = temp_out[j + Ny * i];
		}}
	}
	free(temp_in); free(temp_out);
}




























	
// Mulithreaded Version --> FASTER :)
typedef struct
{
	
	int Nx; 
	int Ny;
	int Nz;
	double complex * in_array;
	double complex * out_array;
	int NUM_THREADS;
	
	int thread_num;
	
} IFFT2_slices_MultiThreadedStruct;

// Helper function that each thread runs
// Goes from f * i -> f * i - 1 in vector calculation
void * IFFT2_slices_MultiThreaded_Helper(void * arg)
{
	IFFT2_slices_MultiThreadedStruct * actual_arg = 
		(IFFT2_slices_MultiThreadedStruct *) arg;
	int Nx = actual_arg -> Nx; 
	int Ny = actual_arg -> Ny;
	int Nz = actual_arg -> Nz;
	double complex * in_array = actual_arg -> in_array;
	double complex * out_array = actual_arg -> out_array;
	int NUM_THREADS = actual_arg -> NUM_THREADS;
	
	int thread_num = actual_arg -> thread_num;

	// Where to calculate?
	int total = Nz;
	int f = (int) total / NUM_THREADS;
	int start = thread_num * f;
	int end = (thread_num + 1) * f;
	if(end >= total)
		end = total;

	IFFT2_slices(Nx, Ny, Nz, in_array, out_array, start, end);	

	return NULL;
}

// Main driver
void IFFT2_slices_MultiThreaded(int Nx, int Ny, int Nz, 
	double complex * in_array, double complex * out_array, int NUM_THREADS)
{
	pthread_t thread[NUM_THREADS];
	int tid[NUM_THREADS];
	IFFT2_slices_MultiThreadedStruct args[NUM_THREADS];
	
	for(int i = 0; i < NUM_THREADS; i++)
	{
		tid[i] = i;
		
		// argument construction
		args[i].Nx = Nx;
		args[i].Ny = Ny;
		args[i].Nz = Nz;

		args[i].in_array = in_array;
		args[i].out_array = out_array;
		args[i].NUM_THREADS = NUM_THREADS;
		
		args[i].thread_num = i;

		// create thread
		pthread_create(&thread[i], NULL, IFFT2_slices_MultiThreaded_Helper, &args[i]);
	}
	
	// wait for threads to finish
	for(int i = 0; i < NUM_THREADS; i++)
		pthread_join(thread[i], NULL);
}

/***
*
* --------------------------------------------------------------------------
*
***/