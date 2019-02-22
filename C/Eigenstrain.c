#include "Eigenstrain.h"

void InitEigenstrain(int Nx, int Ny, int Nz, struct EigenstrainStruct * Eigenstrain)
{
	Eigenstrain -> Eigenstrain_11_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_22_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_33_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_12_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));	
	Eigenstrain -> Eigenstrain_13_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_23_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	Eigenstrain -> Eigenstrain_11_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_22_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_33_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_12_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));	
	Eigenstrain -> Eigenstrain_13_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	Eigenstrain -> Eigenstrain_23_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));	
}


void CalcEigenstrain(struct EigenstrainStruct * Eigenstrain, 
	struct PVector * P, struct ConstantsStruct * Constants)
{
	for(int i = 0; i < Constants -> total; i++)
	{
		Eigenstrain -> Eigenstrain_11_RowMaj[i] = 
			Constants->Q11 * P->P1_RowMaj[i] * P->P1_RowMaj[i] + 
			Constants->Q12 * (P->P2_RowMaj[i] * P->P2_RowMaj[i] + 
							  P->P3_RowMaj[i] * P->P3_RowMaj[i]);

		Eigenstrain -> Eigenstrain_22_RowMaj[i] = 
			Constants->Q11 * P->P2_RowMaj[i] * P->P2_RowMaj[i] + 
			Constants->Q12 * (P->P1_RowMaj[i] * P->P1_RowMaj[i] + 
						      P->P3_RowMaj[i] * P->P3_RowMaj[i]);

		Eigenstrain -> Eigenstrain_33_RowMaj[i] = 
			Constants->Q11 * P->P3_RowMaj[i] * P->P3_RowMaj[i] + 
			Constants->Q12 * (P->P1_RowMaj[i] * P->P1_RowMaj[i] + 
							  P->P2_RowMaj[i] * P->P2_RowMaj[i]);
		
		Eigenstrain -> Eigenstrain_23_RowMaj[i] = 
			Constants->Q44 * P->P2_RowMaj[i] * P->P3_RowMaj[i];
		
		Eigenstrain -> Eigenstrain_13_RowMaj[i] = 
			Constants->Q44 * P->P1_RowMaj[i] * P->P3_RowMaj[i];
		
		Eigenstrain -> Eigenstrain_12_RowMaj[i] = 
			Constants->Q44 * P->P1_RowMaj[i] * P->P2_RowMaj[i];
	}
	
	FFT3(Constants->Nx, Constants->Ny, Constants->Nz, 
		Eigenstrain -> Eigenstrain_11_RowMaj, Eigenstrain -> Eigenstrain_11_3Dk_RowMaj);
	
	FFT3(Constants->Nx, Constants->Ny, Constants->Nz, 
		Eigenstrain -> Eigenstrain_22_RowMaj, Eigenstrain -> Eigenstrain_22_3Dk_RowMaj);
	
	FFT3(Constants->Nx, Constants->Ny, Constants->Nz, 
		Eigenstrain -> Eigenstrain_33_RowMaj, Eigenstrain -> Eigenstrain_33_3Dk_RowMaj);
	
	FFT3(Constants->Nx, Constants->Ny, Constants->Nz, 
		Eigenstrain -> Eigenstrain_12_RowMaj, Eigenstrain -> Eigenstrain_12_3Dk_RowMaj);
	
	FFT3(Constants->Nx, Constants->Ny, Constants->Nz, 
		Eigenstrain -> Eigenstrain_13_RowMaj, Eigenstrain -> Eigenstrain_13_3Dk_RowMaj);
	
	FFT3(Constants->Nx, Constants->Ny, Constants->Nz, 
		Eigenstrain -> Eigenstrain_23_RowMaj, Eigenstrain -> Eigenstrain_23_3Dk_RowMaj);

}


void FreeEigenstrain(struct EigenstrainStruct * Eigenstrain)
{
	free(Eigenstrain -> Eigenstrain_11_RowMaj);
	free(Eigenstrain -> Eigenstrain_22_RowMaj);
	free(Eigenstrain -> Eigenstrain_33_RowMaj);
	free(Eigenstrain -> Eigenstrain_12_RowMaj);	
	free(Eigenstrain -> Eigenstrain_13_RowMaj);
	free(Eigenstrain -> Eigenstrain_23_RowMaj);
	
	free(Eigenstrain -> Eigenstrain_11_3Dk_RowMaj);
	free(Eigenstrain -> Eigenstrain_22_3Dk_RowMaj);
	free(Eigenstrain -> Eigenstrain_33_3Dk_RowMaj);
	free(Eigenstrain -> Eigenstrain_12_3Dk_RowMaj);	
	free(Eigenstrain -> Eigenstrain_13_3Dk_RowMaj);
	free(Eigenstrain -> Eigenstrain_23_3Dk_RowMaj);		
}