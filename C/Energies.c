#include "Energies.h"

void FreeEnergies(struct EnergiesStruct * F)
{

	free(F -> f1_RowMaj);
	free(F -> f2_RowMaj);
	free(F -> f3_RowMaj);

	free(F -> f1_3Dk_RowMaj);
	free(F -> f2_3Dk_RowMaj);
	free(F -> f3_3Dk_RowMaj);	
}

void InitEnergies(int Nx, int Ny, int Nz, struct EnergiesStruct * F)
{
	F -> f1_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	F -> f2_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	F -> f3_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

	F -> f1_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	F -> f2_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	F -> f3_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
}