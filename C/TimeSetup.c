#include "TimeSetup.h"

void TimeSetup( struct ConstantsStruct * Constants, 
				struct TimeSolutionStruct * TimeSolution)
{
	// Load Constants
	double dt = Constants -> dt;
	int Nx = Constants -> Nx;
	int Ny = Constants -> Ny;
	double h_int = Constants -> h_int;
	double h_film = Constants -> h_film;
	double G11 = Constants -> G11;
	double H44 = Constants -> H44;
	double ** kx_grid_2D = Constants -> kx_grid_2D;
	double ** ky_grid_2D = Constants -> ky_grid_2D;
	
	// First Order Calculations
	double complex ** time_q_1_1storder = AllocLargeArray_2D_Complex(Nx, Ny);
	double complex ** time_q_2_1storder = AllocLargeArray_2D_Complex(Nx, Ny);
	double complex ** time_q_3_1storder = AllocLargeArray_2D_Complex(Nx, Ny);
	
	double complex **** time_bc_mat_inv_1_1storder = AllocLargeArray_4D_Complex(Nx, Ny, 2, 2);
	double complex **** time_bc_mat_inv_2_1storder = AllocLargeArray_4D_Complex(Nx, Ny, 2, 2);
	double complex **** time_bc_mat_inv_3_1storder = AllocLargeArray_4D_Complex(Nx, Ny, 2, 2);
	
	// Second Order Calculations
	double complex ** time_q_1_2ndorder = AllocLargeArray_2D_Complex(Nx, Ny);
	double complex ** time_q_2_2ndorder = AllocLargeArray_2D_Complex(Nx, Ny);
	double complex ** time_q_3_2ndorder = AllocLargeArray_2D_Complex(Nx, Ny);
		
	double complex **** time_bc_mat_inv_1_2ndorder = AllocLargeArray_4D_Complex(Nx, Ny, 2, 2);
	double complex **** time_bc_mat_inv_2_2ndorder = AllocLargeArray_4D_Complex(Nx, Ny, 2, 2);
	double complex **** time_bc_mat_inv_3_2ndorder = AllocLargeArray_4D_Complex(Nx, Ny, 2, 2);
		
	double eta_1, eta_2;
	double q_1, q_2, q_3;
	double factor;
	
	// First Order Setups
	// BC at h_int and h_film
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
		
		eta_1 = kx_grid_2D[k1][k2];
		eta_2 = ky_grid_2D[k1][k2];
		
		q_1 = sqrt( (1.0 + dt * H44 * eta_2 * eta_2 + dt * G11 * eta_1 * eta_1) / (dt * H44) );
		q_2 = sqrt( (1.0 + dt * H44 * eta_1 * eta_1 + dt * G11 * eta_2 * eta_2) / (dt * H44) );
		q_3 = sqrt( (1.0 + dt * H44 * eta_1 * eta_2 + dt * H44 * eta_2 * eta_2) / (dt * G11) );
		
		time_q_1_1storder[k1][k2] = q_1;
		time_q_2_1storder[k1][k2] = q_2;
		time_q_3_1storder[k1][k2] = q_3;
		
		factor = 1 / (exp(q_1 * h_int) * exp(-q_1 * h_film) - 
					  exp(-q_1 * h_int) * exp(q_1 * h_film) );
		time_bc_mat_inv_1_1storder[k1][k2][0][0] = factor * exp(-q_1 * h_film);   time_bc_mat_inv_1_1storder[k1][k2][0][1] = -factor * exp(-q_1 * h_int);
		time_bc_mat_inv_1_1storder[k1][k2][1][0] = -factor * exp(q_1 * h_film);   time_bc_mat_inv_1_1storder[k1][k2][1][1] = factor * exp(q_1 * h_int);

		factor = 1 / (exp(q_2 * h_int) * exp(-q_2 * h_film) - 
					  exp(-q_2 * h_int) * exp(q_2 * h_film) );		
		time_bc_mat_inv_2_1storder[k1][k2][0][0] = factor * exp(-q_2 * h_film);   time_bc_mat_inv_2_1storder[k1][k2][0][1] = -factor * exp(-q_2 * h_int);
		time_bc_mat_inv_2_1storder[k1][k2][1][0] = -factor * exp(q_2 * h_film);   time_bc_mat_inv_2_1storder[k1][k2][1][1] = factor * exp(q_2 * h_int);
				
		factor = 1 / (exp(q_3 * h_int) * exp(-q_3 * h_film) - 
					  exp(-q_3 * h_int) * exp(q_3 * h_film) );		
		time_bc_mat_inv_3_1storder[k1][k2][0][0] = factor * exp(-q_3 * h_film);   time_bc_mat_inv_3_1storder[k1][k2][0][1] = -factor * exp(-q_3 * h_int);
		time_bc_mat_inv_3_1storder[k1][k2][1][0] = -factor * exp(q_3 * h_film);   time_bc_mat_inv_3_1storder[k1][k2][1][1] = factor * exp(q_3 * h_int);
	}}
	
	// Second order setups
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
		
		eta_1 = kx_grid_2D[k1][k2];
		eta_2 = ky_grid_2D[k1][k2];
		
		q_1 = sqrt( (3 + 2 * dt * G11 * eta_1 * eta_1 + 2 * dt * H44 * eta_2 * eta_2) / 
				(2 * dt * H44) );
		q_1 = sqrt( (3 + 2 * dt * G11 * eta_2 * eta_2 + 2 * dt * H44 * eta_1 * eta_1) / 
				(2 * dt * H44) );
		q_1 = sqrt( (3 + 2 * dt * H44 * eta_1 * eta_1 + 2 * dt * H44 * eta_2 * eta_2) / 
				(2 * dt * G11) );
				
		time_q_1_2ndorder[k1][k2] = q_1;
		time_q_2_2ndorder[k1][k2] = q_2;
		time_q_3_2ndorder[k1][k2] = q_3;
		
		factor = 1 / (exp(q_1 * h_int) * exp(-q_1 * h_film) - 
					  exp(-q_1 * h_int) * exp(q_1 * h_film) );
		time_bc_mat_inv_1_2ndorder[k1][k2][0][0] = factor * exp(-q_1 * h_film);   time_bc_mat_inv_1_2ndorder[k1][k2][0][1] = -factor * exp(-q_1 * h_int);
		time_bc_mat_inv_1_2ndorder[k1][k2][1][0] = -factor * exp(q_1 * h_film);   time_bc_mat_inv_1_2ndorder[k1][k2][1][1] = factor * exp(q_1 * h_int);

		factor = 1 / (exp(q_2 * h_int) * exp(-q_2 * h_film) - 
					  exp(-q_2 * h_int) * exp(q_2 * h_film) );		
		time_bc_mat_inv_2_2ndorder[k1][k2][0][0] = factor * exp(-q_2 * h_film);   time_bc_mat_inv_2_2ndorder[k1][k2][0][1] = -factor * exp(-q_2 * h_int);
		time_bc_mat_inv_2_2ndorder[k1][k2][1][0] = -factor * exp(q_2 * h_film);   time_bc_mat_inv_2_2ndorder[k1][k2][1][1] = factor * exp(q_2 * h_int);
				
		factor = 1 / (exp(q_3 * h_int) * exp(-q_3 * h_film) - 
					  exp(-q_3 * h_int) * exp(q_3 * h_film) );		
		time_bc_mat_inv_3_2ndorder[k1][k2][0][0] = factor * exp(-q_3 * h_film);   time_bc_mat_inv_3_2ndorder[k1][k2][0][1] = -factor * exp(-q_3 * h_int);
		time_bc_mat_inv_3_2ndorder[k1][k2][1][0] = -factor * exp(q_3 * h_film);   time_bc_mat_inv_3_2ndorder[k1][k2][1][1] = factor * exp(q_3 * h_int);
	}}
	
	// Convert...
	//		time_bc_mat_inv_1_1storder
	//		time_bc_mat_inv_2_1storder
	//		time_bc_mat_inv_3_1storder
	//		time_bc_mat_inv_1_2ndorder
	//		time_bc_mat_inv_2_2ndorder
	//		time_bc_mat_inv_3_2ndorder
	double complex *** time_bc_mat_inv_1_1storder_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 4);
	double complex *** time_bc_mat_inv_2_1storder_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 4);
	double complex *** time_bc_mat_inv_3_1storder_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 4);	

	double complex *** time_bc_mat_inv_1_2ndorder_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 4);
	double complex *** time_bc_mat_inv_2_2ndorder_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 4);
	double complex *** time_bc_mat_inv_3_2ndorder_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 4);
	
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		Convert2D_2_ColMajor_Complex(2, 2, 
			time_bc_mat_inv_1_1storder[i][j], 
			time_bc_mat_inv_1_1storder_ColMaj[i][j]);
		Convert2D_2_ColMajor_Complex(2, 2, 
			time_bc_mat_inv_2_1storder[i][j], 
			time_bc_mat_inv_2_1storder_ColMaj[i][j]);
		Convert2D_2_ColMajor_Complex(2, 2, 
			time_bc_mat_inv_3_1storder[i][j], 
			time_bc_mat_inv_3_1storder_ColMaj[i][j]);
			
		Convert2D_2_ColMajor_Complex(2, 2, 
			time_bc_mat_inv_1_2ndorder[i][j], 
			time_bc_mat_inv_1_2ndorder_ColMaj[i][j]);	
		Convert2D_2_ColMajor_Complex(2, 2, 
			time_bc_mat_inv_2_2ndorder[i][j], 
			time_bc_mat_inv_2_2ndorder_ColMaj[i][j]);
		Convert2D_2_ColMajor_Complex(2, 2, 
			time_bc_mat_inv_3_2ndorder[i][j], 
			time_bc_mat_inv_3_2ndorder_ColMaj[i][j]);
	}}
	
	TimeSolution -> time_q_1_1storder = time_q_1_1storder;
	TimeSolution -> time_q_2_1storder = time_q_2_1storder;
	TimeSolution -> time_q_3_1storder = time_q_3_1storder;
	
	TimeSolution -> time_bc_mat_inv_1_1storder_ColMaj = time_bc_mat_inv_1_1storder_ColMaj;
	TimeSolution -> time_bc_mat_inv_2_1storder_ColMaj = time_bc_mat_inv_2_1storder_ColMaj;
	TimeSolution -> time_bc_mat_inv_3_1storder_ColMaj = time_bc_mat_inv_3_1storder_ColMaj;
	
	TimeSolution -> time_q_1_2ndorder = time_q_1_2ndorder;
	TimeSolution -> time_q_2_2ndorder = time_q_2_2ndorder;
	TimeSolution -> time_q_3_2ndorder = time_q_3_2ndorder;
	
	TimeSolution -> time_bc_mat_inv_1_2ndorder_ColMaj = time_bc_mat_inv_1_2ndorder_ColMaj;
	TimeSolution -> time_bc_mat_inv_2_2ndorder_ColMaj = time_bc_mat_inv_2_2ndorder_ColMaj;
	TimeSolution -> time_bc_mat_inv_3_2ndorder_ColMaj = time_bc_mat_inv_3_2ndorder_ColMaj;

	// Free
	Free_4D_mat_Complex(Nx, Ny, 2, 2, time_bc_mat_inv_1_1storder);
	Free_4D_mat_Complex(Nx, Ny, 2, 2, time_bc_mat_inv_2_1storder);
	Free_4D_mat_Complex(Nx, Ny, 2, 2, time_bc_mat_inv_3_1storder);
	
	Free_4D_mat_Complex(Nx, Ny, 2, 2, time_bc_mat_inv_1_2ndorder);
	Free_4D_mat_Complex(Nx, Ny, 2, 2, time_bc_mat_inv_2_2ndorder);
	Free_4D_mat_Complex(Nx, Ny, 2, 2, time_bc_mat_inv_3_2ndorder);	
	
}


void FreeTimeSolution( int Nx, int Ny, struct TimeSolutionStruct * TimeSolution)
{
	Free_2D_mat_Complex(Nx, Ny, TimeSolution -> time_q_1_1storder);
	Free_2D_mat_Complex(Nx, Ny, TimeSolution -> time_q_2_1storder);
	Free_2D_mat_Complex(Nx, Ny, TimeSolution -> time_q_3_1storder);
	
	Free_3D_mat_Complex(Nx, Ny, 4, TimeSolution -> time_bc_mat_inv_1_1storder_ColMaj);
	Free_3D_mat_Complex(Nx, Ny, 4, TimeSolution -> time_bc_mat_inv_2_1storder_ColMaj);
	Free_3D_mat_Complex(Nx, Ny, 4, TimeSolution -> time_bc_mat_inv_3_1storder_ColMaj);

	Free_2D_mat_Complex(Nx, Ny, TimeSolution -> time_q_1_2ndorder);
	Free_2D_mat_Complex(Nx, Ny, TimeSolution -> time_q_2_2ndorder);
	Free_2D_mat_Complex(Nx, Ny, TimeSolution -> time_q_3_2ndorder);
	
	Free_3D_mat_Complex(Nx, Ny, 4, TimeSolution -> time_bc_mat_inv_1_2ndorder_ColMaj);
	Free_3D_mat_Complex(Nx, Ny, 4, TimeSolution -> time_bc_mat_inv_2_2ndorder_ColMaj);
	Free_3D_mat_Complex(Nx, Ny, 4, TimeSolution -> time_bc_mat_inv_3_2ndorder_ColMaj);

}