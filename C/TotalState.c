#include "TotalState.h"

void FreeTotalState(struct TotalStateStruct * TotalState)
{
	free(TotalState -> Potential_RowMaj);

	free(TotalState -> E_1_depol_RowMaj);
	free(TotalState -> E_2_depol_RowMaj);
	free(TotalState -> E_3_depol_RowMaj);

	free(TotalState -> f1_elec_RowMaj);
	free(TotalState -> f2_elec_RowMaj);
	free(TotalState -> f3_elec_RowMaj);	
	
	free(TotalState -> f1_elastic_RowMaj);
	free(TotalState -> f2_elastic_RowMaj);
	free(TotalState -> f3_elastic_RowMaj);
	
	free(TotalState -> u_1_RowMaj);
	free(TotalState -> u_2_RowMaj);
	free(TotalState -> u_3_RowMaj);
	
	free(TotalState -> e_11_RowMaj);
	free(TotalState -> e_12_RowMaj);
	free(TotalState -> e_13_RowMaj);
	free(TotalState -> e_22_RowMaj);
	free(TotalState -> e_23_RowMaj);
	free(TotalState -> e_33_RowMaj);
	
	free(TotalState -> TotalStrain_11_RowMaj);
	free(TotalState -> TotalStrain_12_RowMaj);
	free(TotalState -> TotalStrain_13_RowMaj);
	free(TotalState -> TotalStrain_22_RowMaj);
	free(TotalState -> TotalStrain_23_RowMaj);
	free(TotalState -> TotalStrain_33_RowMaj);
	
	free(TotalState -> f1_landau_RowMaj);
	free(TotalState -> f2_landau_RowMaj);
	free(TotalState -> f3_landau_RowMaj);
	
}

void InitTotalState(int Nx, int Ny, int Nz, struct TotalStateStruct * TotalState)
{
	TotalState -> Potential_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

	TotalState -> E_1_depol_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> E_2_depol_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> E_3_depol_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

	TotalState -> f1_elec_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> f2_elec_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> f3_elec_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	TotalState -> f1_elastic_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> f2_elastic_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> f3_elastic_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	TotalState -> u_1_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> u_2_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> u_3_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	TotalState -> e_11_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> e_12_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> e_13_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> e_22_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> e_23_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> e_33_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	TotalState -> TotalStrain_11_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> TotalStrain_12_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> TotalStrain_13_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> TotalStrain_22_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> TotalStrain_23_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> TotalStrain_33_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	TotalState -> f1_landau_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> f2_landau_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	TotalState -> f3_landau_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
}