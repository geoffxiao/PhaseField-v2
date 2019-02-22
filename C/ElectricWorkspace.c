#include "ElectricWorkspace.h"

void FreeElectricWorkspace(struct ElectricWorkspaceStruct * E)
{
	free(E -> potential_A_3Dk_RowMaj);
	free(E -> potential_A_RowMaj);

	free(E -> potential_B_d3_2Dk_RowMaj);
	free(E -> potential_B_2Dk_RowMaj);
	free(E -> potential_B_RowMaj);
		
	free(E -> E_A_1_3Dk_RowMaj);
	free(E -> E_A_2_3Dk_RowMaj);
	free(E -> E_A_3_3Dk_RowMaj);
	free(E -> E_A_1_RowMaj);
	free(E -> E_A_2_RowMaj);
	free(E -> E_A_3_RowMaj);	
	
	free(E -> E_B_1_2Dk_RowMaj);
	free(E -> E_B_2_2Dk_RowMaj);
	free(E -> E_B_3_2Dk_RowMaj);
	
	free(E -> E_B_1_RowMaj);
	free(E -> E_B_2_RowMaj);
	free(E -> E_B_3_RowMaj);
		
	free(E -> potential_bc_interface_RowMaj);
	free(E -> potential_bc_film_RowMaj);
	
	free(E -> potential_bc_interface_2Dk_RowMaj);
	free(E -> potential_bc_film_2Dk_RowMaj);
}

void InitElectricWorkspace(int Nx, int Ny, int Nz, struct ElectricWorkspaceStruct * E)
{
	E -> potential_A_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E -> potential_A_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	E -> potential_B_2Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E -> potential_B_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E -> potential_B_d3_2Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	E ->  E_B_1_2Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E ->  E_B_2_2Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E ->  E_B_3_2Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	E ->  E_B_1_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E ->  E_B_2_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E ->  E_B_3_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	E -> E_A_1_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	E -> E_A_2_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));	
	E -> E_A_3_3Dk_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));	

	E -> E_A_1_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));	
	E -> E_A_2_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));	
	E -> E_A_3_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));

	E -> potential_bc_interface_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	E -> potential_bc_film_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	
	E -> potential_bc_interface_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	E -> potential_bc_film_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));


}