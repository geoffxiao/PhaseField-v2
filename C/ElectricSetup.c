#include "ElectricSetup.h"

void ElectricSetup( struct ConstantsStruct * Constants,
					struct ElectricSolutionStruct * ElectricSolution)
{
	// Extract Constants from ConstantsStruct
	int Nx = Constants -> Nx;
	int Ny = Constants -> Ny;
	double h_film = Constants -> h_film;
	double h_int = Constants -> h_int;
	double k_electric_11 = Constants -> k_electric_11;
	double k_electric_22 = Constants -> k_electric_22;
	double k_electric_33 = Constants -> k_electric_33;
	double ** kx_grid_2D = Constants -> kx_grid_2D;
	double ** ky_grid_2D = Constants -> ky_grid_2D;
	
	// Allocate Memory for the Outputs
	double **** electric_bc_mat = AllocLargeArray_4D(Nx, Ny, 2, 2);
	double **** electric_bc_mat_inv = AllocLargeArray_4D(Nx, Ny, 2, 2);
	double ** electric_p_mat = AllocLargeArray_2D(Nx, Ny);
	double ** electric_bc_mat_inv_korigin = AllocLargeArray_2D(2, 2);
	double ** potential_sol_mat = AllocLargeArray_2D(2, 2);
	
	double eta_1, eta_2, p, factor;
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
	if(k1 != 0 || k2 != 0){
		
		eta_1 = kx_grid_2D[k1][k2];
		eta_2 = ky_grid_2D[k1][k2];
		
		p = sqrt( (k_electric_11 * eta_1 * eta_1 +
				   k_electric_22 * eta_2 * eta_2) / k_electric_33 );

		potential_sol_mat[0][0] = exp(h_int * p);    potential_sol_mat[0][1] = exp(-h_int * p);
		potential_sol_mat[1][0] = exp(h_film * p);   potential_sol_mat[1][1] = exp(-h_film * p);
		
		electric_bc_mat[k1][k2][0][0] = potential_sol_mat[0][0];   electric_bc_mat[k1][k2][0][1] = potential_sol_mat[0][1];
		electric_bc_mat[k1][k2][1][0] = potential_sol_mat[1][0];   electric_bc_mat[k1][k2][1][1] = potential_sol_mat[1][1];

		factor = ( 1 / (potential_sol_mat[0][0] * potential_sol_mat[1][1] - 
						potential_sol_mat[0][1] * potential_sol_mat[1][0]) );
		electric_bc_mat_inv[k1][k2][0][0] = factor * potential_sol_mat[1][1];   electric_bc_mat_inv[k1][k2][0][1] = -factor * potential_sol_mat[0][1];
		electric_bc_mat_inv[k1][k2][1][0] = -factor * potential_sol_mat[1][0];  electric_bc_mat_inv[k1][k2][1][1] = factor * potential_sol_mat[0][0];
		
		electric_p_mat[k1][k2] = p;
	}}}
	
	factor = 1 / (h_int - h_film);
	electric_bc_mat_inv_korigin[0][0] = factor;             electric_bc_mat_inv_korigin[0][1] = -factor;
	electric_bc_mat_inv_korigin[1][0] = -factor * h_film;   electric_bc_mat_inv_korigin[1][1] = factor * h_int;

	
	// Convert...
	// 		electric_bc_mat_inv
	//		electric_bc_mat_inv_korigin
	double complex *** electric_bc_mat_inv_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 4);
	double complex * electric_bc_mat_inv_korigin_ColMaj = malloc(sizeof(double complex) * 4);
	
	Convert2D_2_ColMajor_Double_2_Complex(2, 2, 
		electric_bc_mat_inv_korigin, 
		electric_bc_mat_inv_korigin_ColMaj);
		
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		Convert2D_2_ColMajor_Double_2_Complex(2, 2, 
			electric_bc_mat_inv[i][j], 
			electric_bc_mat_inv_ColMaj[i][j]);
	}}	
	// Free memory
	Free_2D_mat(2, 2, potential_sol_mat);
	Free_2D_mat(2, 2, electric_bc_mat_inv_korigin);
	Free_4D_mat(Nx, Ny, 2, 2, electric_bc_mat_inv);
	Free_4D_mat(Nx, Ny, 2, 2, electric_bc_mat);
	
	// Outputs
	ElectricSolution -> electric_bc_mat_inv_korigin_ColMaj = electric_bc_mat_inv_korigin_ColMaj;
	ElectricSolution -> electric_bc_mat_inv_ColMaj = electric_bc_mat_inv_ColMaj;
	ElectricSolution -> electric_p_mat = electric_p_mat;
	// ElectricSolution -> electric_bc_mat = electric_bc_mat;

}


void FreeElectricSolution( int Nx, int Ny, struct ElectricSolutionStruct * ElectricSolution)
{
	Free_3D_mat_Complex(Nx, Ny, 4, ElectricSolution -> electric_bc_mat_inv_ColMaj);
	free(ElectricSolution -> electric_bc_mat_inv_korigin_ColMaj);
	Free_2D_mat(Nx, Ny, ElectricSolution -> electric_p_mat);
	// Free_4D_mat(Nx, Ny, 2, 2, ElectricSolution -> electric_bc_mat);
}