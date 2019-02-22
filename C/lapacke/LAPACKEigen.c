#include <stdio.h>
#include <math.h>

#include <complex.h>
#include <stdlib.h>

//.......................................................................................................
void zgeTranspose( double complex *Transposed, double complex *M ,int n)
{
	int i,j;
	for(i=0;i<n;i++)
	for(j=0;j<n;j++) 
		Transposed[i+n*j] = M[i*n+j];
}

//......................................................................................................
//  MatrixComplexEigensystem: computes the eigenvectors and eigenValues of input matrix A
//  The eigenvectors are stored in columns
//.....................................................................................................
void MatrixComplexEigensystem( double complex *eigenvectorsVR, double complex *eigenvaluesW, double complex *A, int N)
{
	int i;
	double complex *AT = (double complex*) malloc( N*N*sizeof(double complex) );
	zgeTranspose( AT, A , N);
	char JOBVL ='N';   // Compute Right eigenvectors
	char JOBVR ='V';   // Do not compute Left eigenvectors
	double complex VL[1];

	int LDVL = 1; 
	int LDVR = N;
	int LWORK = 4*N; 

	double complex *WORK =  (double complex*)malloc( LWORK*sizeof(double complex));
	double complex *RWORK = (double complex*)malloc( 2*N*sizeof(double complex));

	int INFO;

	zgeev_( &JOBVL, &JOBVR, &N, AT ,  &N , eigenvaluesW ,

	VL, &LDVL,
	eigenvectorsVR, &LDVR, 
	WORK, 
	&LWORK, RWORK, &INFO );

	zgeTranspose( AT, eigenvectorsVR , N);

	for(i=0;i<N*N;i++) 
		eigenvectorsVR[i]=AT[i];
	free(WORK);
	free(RWORK);
	free(AT);
}


int main()
{
	int i,j;
	const int N = 3;

	double complex A[] = { 1.+I , 2. ,  3 , 4. , 5.+I , 6. , 7., 8., 9. + I};

	double complex eigenVectors[N*N];
	double complex eigenValues[N];


	MatrixComplexEigensystem( eigenVectors, eigenValues, A, N);

	printf("\nEigenvectors\n");

	for(i=0;i<N;i++){
		for(j=0;j<N;j++) 
			printf(" (%f,%f) \t", eigenVectors[i*N + j]);
		printf("\n");
	}

	printf("\nEigenvalues \n");
	for(i=0;i<N;i++) 
		printf("\n (%f, %f) \t",  eigenValues[i] );


	printf("\n------------------------------------------------------------\n"); 
	return 0;
}