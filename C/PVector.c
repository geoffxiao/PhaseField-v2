#include "PVector.h"

void FreePVector(struct PVector * P)
{

	free(P -> P1_RowMaj);
	free(P -> P2_RowMaj);
	free(P -> P3_RowMaj);
	
	free(P -> P1_3Dk_RowMaj);
	free(P -> P2_3Dk_RowMaj);
	free(P -> P3_3Dk_RowMaj);	
}

// type = -1 to zero P
void InitPVector(int Nx, int Ny, int Nz, int type, struct PVector * P)
{
	if(type == 1) // load from file
	{
		P -> P1_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P2_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P3_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

		FILE * P1_fd = fopen("P1.txt", "r");
		FILE * P2_fd = fopen("P2.txt", "r");
		FILE * P3_fd = fopen("P3.txt", "r");
		
		double a, b, c;
		// Randomize initial conditions
		for(int i = 0; i < Nx * Ny * Nz; i++)
		{
			fscanf(P1_fd, "%lf", &a);
			fscanf(P2_fd, "%lf", &b);
			fscanf(P3_fd, "%lf", &c);

			P -> P1_RowMaj[i] = a;
			P -> P2_RowMaj[i] = b;
			P -> P3_RowMaj[i] = c;
		}
		fclose(P1_fd); fclose(P2_fd); fclose(P3_fd);
		
		P -> P1_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P2_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P3_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

		FFT3(Nx, Ny, Nz, P -> P1_RowMaj, P -> P1_3Dk_RowMaj);
		FFT3(Nx, Ny, Nz, P -> P2_RowMaj, P -> P2_3Dk_RowMaj);
		FFT3(Nx, Ny, Nz, P -> P3_RowMaj, P -> P3_3Dk_RowMaj);
	}
	else if(type == 0) // don't load from file
	{
		P -> P1_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P2_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P3_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

		for(int i = 0; i < Nx * Ny * Nz; i++)
		{
			P -> P1_RowMaj[i] = ((double)rand()/RAND_MAX*2.0-1.0) * 0.001;
			P -> P2_RowMaj[i] = ((double)rand()/RAND_MAX*2.0-1.0) * 0.001;
			P -> P3_RowMaj[i] = ((double)rand()/RAND_MAX*2.0-1.0) * 0.001;
		}		
		P -> P1_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P2_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P3_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

		FFT3(Nx, Ny, Nz, P -> P1_RowMaj, P -> P1_3Dk_RowMaj);
		FFT3(Nx, Ny, Nz, P -> P2_RowMaj, P -> P2_3Dk_RowMaj);
		FFT3(Nx, Ny, Nz, P -> P3_RowMaj, P -> P3_3Dk_RowMaj);
	}
	else if(type == 2)
	{
		P -> P1_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P2_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P3_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

		P -> P1_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P2_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		P -> P3_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		
	}
}