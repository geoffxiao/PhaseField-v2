#include "CalcElasticEnergy.h"

void CalcElasticEnergy(struct ElasticWorkspace * U, 
	struct PVector * P,
	struct EigenstrainStruct * ES,
	struct ConstantsStruct * C, 
	struct ElasticSolutionStruct * S,
	struct TotalStateStruct * T)
{
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	CalcEigenstrain(ES, P, C);
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	
	double complex * u_1_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_2_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_3_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));

	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	// Particular solution
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	for(int i = 0; i < C->total; i++)
	{
		u_1_A_3Dk_RowMaj[i] = 
		  - C->C11*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_11_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_11_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_11_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_12_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_12_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_13_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_13_RowMaj[i]*2i - 
			C->C12*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_12_RowMaj[i]*I - 
			C->C11*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_12_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_12_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_12_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_11_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_23_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_13_RowMaj[i]*2i - 
			C->C12*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_13_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_13_RowMaj[i]*I - 
			C->C11*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_13_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_13_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_11_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_23_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_12_RowMaj[i]*2i;
		u_2_A_3Dk_RowMaj[i] = 
		  - C->C11*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_21_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_21_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_21_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_12_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_22_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_13_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_23_RowMaj[i]*2i - 
			C->C12*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_22_RowMaj[i]*I - 
			C->C11*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_22_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_22_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_12_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_21_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_23_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_23_RowMaj[i]*2i - 
			C->C12*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_23_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_23_RowMaj[i]*I - 
			C->C11*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_23_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_13_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_21_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_23_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_22_RowMaj[i]*2i;
		u_3_A_3Dk_RowMaj[i] = 
		  - C->C11*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_31_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_31_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_31_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_12_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_32_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_13_3Dk_RowMaj[i]*C->kx_grid_3D_RowMaj[i]*S->Green_33_RowMaj[i]*2i - 
			C->C12*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_32_RowMaj[i]*I - 
			C->C11*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_32_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_32_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_12_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_31_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_23_3Dk_RowMaj[i]*C->ky_grid_3D_RowMaj[i]*S->Green_33_RowMaj[i]*2i - 
			C->C12*ES->Eigenstrain_11_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_33_RowMaj[i]*I - 
			C->C12*ES->Eigenstrain_22_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_33_RowMaj[i]*I - 
			C->C11*ES->Eigenstrain_33_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_33_RowMaj[i]*I - 
			C->C44*ES->Eigenstrain_13_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_31_RowMaj[i]*2i - 
			C->C44*ES->Eigenstrain_23_3Dk_RowMaj[i]*C->kz_grid_3D_RowMaj[i]*S->Green_32_RowMaj[i]*2i;	
	}
	
	// k1 = k2 = 0 point
	u_1_A_3Dk_RowMaj[0] = 0;
	u_2_A_3Dk_RowMaj[0] = 0;
	u_3_A_3Dk_RowMaj[0] = 0;
	
	
	// Find strain A
	double complex * e_11_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_12_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_13_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_22_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_23_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_33_A_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	
	for(int i = 0; i < C->total; i++)
	{
		e_11_A_3Dk_RowMaj[i] = I * C->kx_grid_3D_RowMaj[i] * u_1_A_3Dk_RowMaj[i];
		e_22_A_3Dk_RowMaj[i] = I * C->ky_grid_3D_RowMaj[i] * u_2_A_3Dk_RowMaj[i];
		e_33_A_3Dk_RowMaj[i] = I * C->kz_grid_3D_RowMaj[i] * u_3_A_3Dk_RowMaj[i];
		e_12_A_3Dk_RowMaj[i] = 0.5 * ( I * C->kx_grid_3D_RowMaj[i] * u_2_A_3Dk_RowMaj[i] + I * C->ky_grid_3D_RowMaj[i] * u_1_A_3Dk_RowMaj[i] );
		e_13_A_3Dk_RowMaj[i] = 0.5 * ( I * C->kx_grid_3D_RowMaj[i] * u_3_A_3Dk_RowMaj[i] + I * C->kz_grid_3D_RowMaj[i] * u_1_A_3Dk_RowMaj[i] );
		e_23_A_3Dk_RowMaj[i] = 0.5 * ( I * C->kz_grid_3D_RowMaj[i] * u_2_A_3Dk_RowMaj[i] + I * C->ky_grid_3D_RowMaj[i] * u_3_A_3Dk_RowMaj[i] );
	}
	
	
	double complex * u_1_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_2_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_3_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_11_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_12_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_13_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_22_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_23_A_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_33_A_RowMaj = calloc(C->total, sizeof(double complex));
	
	// IFFT to real space
	IFFT3(C->Nx, C->Ny, C->Nz, u_1_A_3Dk_RowMaj, u_1_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, u_2_A_3Dk_RowMaj, u_2_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, u_3_A_3Dk_RowMaj, u_3_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, e_11_A_3Dk_RowMaj, e_11_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, e_12_A_3Dk_RowMaj, e_12_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, e_13_A_3Dk_RowMaj, e_13_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, e_22_A_3Dk_RowMaj, e_22_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, e_23_A_3Dk_RowMaj, e_23_A_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, e_33_A_3Dk_RowMaj, e_33_A_RowMaj);
	
	
	
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	// Find the boundary conditions
	// Film and Sub
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	
	// Film BC in 3Dk
	double complex * bc_film_1_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * bc_film_2_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * bc_film_3_3Dk_RowMaj = calloc(C->total, sizeof(double complex));
	
	for(int i = 0; i < C->total; i++)
	{
		bc_film_1_3Dk_RowMaj[i] = C->C44*(ES->Eigenstrain_13_3Dk_RowMaj[i] - C->kx_grid_3D_RowMaj[i]*u_3_A_3Dk_RowMaj[i]*I) + C->C44*(ES->Eigenstrain_13_3Dk_RowMaj[i] - C->kz_grid_3D_RowMaj[i]*u_1_A_3Dk_RowMaj[i]*I);
		bc_film_2_3Dk_RowMaj[i] = C->C44*(ES->Eigenstrain_23_3Dk_RowMaj[i] - C->ky_grid_3D_RowMaj[i]*u_3_A_3Dk_RowMaj[i]*I) + C->C44*(ES->Eigenstrain_23_3Dk_RowMaj[i] - C->kz_grid_3D_RowMaj[i]*u_2_A_3Dk_RowMaj[i]*I);
		bc_film_3_3Dk_RowMaj[i] = C->C12*(ES->Eigenstrain_11_3Dk_RowMaj[i] - C->kx_grid_3D_RowMaj[i]*u_1_A_3Dk_RowMaj[i]*I) + C->C12*(ES->Eigenstrain_22_3Dk_RowMaj[i] - C->ky_grid_3D_RowMaj[i]*u_2_A_3Dk_RowMaj[i]*I) + C->C11*(ES->Eigenstrain_33_3Dk_RowMaj[i] - C->kz_grid_3D_RowMaj[i]*u_3_A_3Dk_RowMaj[i]*I);
	}	
	
	// Film BC in 3D Real
	double complex * bc_film_1_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * bc_film_2_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * bc_film_3_RowMaj = calloc(C->total, sizeof(double complex));

	IFFT3(C->Nx, C->Ny, C->Nz, bc_film_1_3Dk_RowMaj, bc_film_1_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, bc_film_2_3Dk_RowMaj, bc_film_2_RowMaj);
	IFFT3(C->Nx, C->Ny, C->Nz, bc_film_3_3Dk_RowMaj, bc_film_3_RowMaj);
	
	// Get the film and substrate bc 2D space, one layer of 3D data
	double complex * bc_film_1_2D_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_film_2_2D_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_film_3_2D_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));

	double complex * bc_sub_1_2D_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_sub_2_2D_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_sub_3_2D_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));

	double complex * bc_film_1_2Dk_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_film_2_2Dk_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_film_3_2Dk_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));

	double complex * bc_sub_1_2Dk_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_sub_2_2Dk_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
	double complex * bc_sub_3_2Dk_RowMaj = calloc(C->Nx * C->Ny, sizeof(double complex));
		
	// extract the z = film_index and z = sub_index points in the 2D plane
	int index_RowMajor_3D = 0;
	int index_RowMajor_2D = 0;
	for(int i = 0; i < C->Nx; i++){
	for(int j = 0; j < C->Ny; j++){
		
		index_RowMajor_3D = C->film_index + C->Nz * (j + C->Ny * i);
		index_RowMajor_2D = j + C->Ny * i;
		
		bc_film_1_2D_RowMaj[index_RowMajor_2D] = bc_film_1_RowMaj[index_RowMajor_3D];
		bc_film_2_2D_RowMaj[index_RowMajor_2D] = bc_film_2_RowMaj[index_RowMajor_3D];
		bc_film_3_2D_RowMaj[index_RowMajor_2D] = bc_film_3_RowMaj[index_RowMajor_3D];
		
		index_RowMajor_3D = C->Nz * (j + C->Ny * i);		
		bc_sub_1_2D_RowMaj[index_RowMajor_2D] = -u_1_A_RowMaj[index_RowMajor_3D];
		bc_sub_2_2D_RowMaj[index_RowMajor_2D] = -u_2_A_RowMaj[index_RowMajor_3D];
		bc_sub_3_2D_RowMaj[index_RowMajor_2D] = -u_3_A_RowMaj[index_RowMajor_3D];
		
	}}
	
	// 2D real space BC to 2Dk BC
	FFT2(C->Nx, C->Ny, bc_film_1_2D_RowMaj, bc_film_1_2Dk_RowMaj);	
	FFT2(C->Nx, C->Ny, bc_film_2_2D_RowMaj, bc_film_2_2Dk_RowMaj);	
	FFT2(C->Nx, C->Ny, bc_film_3_2D_RowMaj, bc_film_3_2Dk_RowMaj);	

	FFT2(C->Nx, C->Ny, bc_sub_1_2D_RowMaj, bc_sub_1_2Dk_RowMaj);	
	FFT2(C->Nx, C->Ny, bc_sub_2_2D_RowMaj, bc_sub_2_2Dk_RowMaj);	
	FFT2(C->Nx, C->Ny, bc_sub_3_2D_RowMaj, bc_sub_3_2Dk_RowMaj);	
	
	// rescale the film bc
	for(int i = 0; i < C->Nx * C->Ny; i++)
	{
		double rescale = C->kx_grid_2D_RowMaj[i] * C->kx_grid_2D_RowMaj[i] + 
				 C->ky_grid_2D_RowMaj[i] * C->ky_grid_2D_RowMaj[i];
		if(i == 0)
			rescale = 1;
		
		bc_film_1_2Dk_RowMaj[i] = bc_film_1_2Dk_RowMaj[i] / rescale;
		bc_film_2_2Dk_RowMaj[i] = bc_film_2_2Dk_RowMaj[i] / rescale;
		bc_film_3_2Dk_RowMaj[i] = bc_film_3_2Dk_RowMaj[i] / rescale;
	}

	
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	// --- homo solution --- // 
	// apply the bc
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	// Construct B
	double complex * u_1_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_2_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_3_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_1_B_d3_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_2_B_d3_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_3_B_d3_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	
	// For LAPACK zgemm
	char N = 'N';
	int n6 = 6;
	int n1 = 1;
	double complex one = 1.0;
	double complex zero = 0.0;
	double complex coeffs[6] = {0,0,0,0,0,0};
	double complex bc[6] = {0,0,0,0,0,0};

	double K_vec;
	for(int k1 = 0; k1 < C->Nx; k1++){
	for(int k2 = 0; k2 < C->Ny; k2++){
		
		K_vec = sqrt(C->kx_grid_2D[k1][k2] * C->kx_grid_2D[k1][k2] + 
					 C->ky_grid_2D[k1][k2] * C->ky_grid_2D[k1][k2]);
		
		bc[0] = bc_sub_1_2Dk_RowMaj[k2 + C->Ny * k1];
		bc[1] = bc_sub_2_2Dk_RowMaj[k2 + C->Ny * k1];
		bc[2] = bc_sub_3_2Dk_RowMaj[k2 + C->Ny * k1];
		bc[3] = bc_film_1_2Dk_RowMaj[k2 + C->Ny * k1];
		bc[4] = bc_film_2_2Dk_RowMaj[k2 + C->Ny * k1];
		bc[5] = bc_film_3_2Dk_RowMaj[k2 + C->Ny * k1];
	
		// A = m x k
		// B = k x n
		// C = m x n
		// Inputs: m, n, k
		//	       m = 6, n = 1, k = 6
		zgemm_(&N, &N, &n6, &n1, &n6, 
				&one, S -> elastic_bc_mat_inv_ColMaj[k1][k2], &n6, 
				bc, &n6, 
				&zero, coeffs, &n6);	

		// P * A * x = L * U * x = b
		// U * x = inv(L) * P * b
		
		// A * x = b
		// x = inv(A) * b
		
		double complex p_1 = S->elastic_eigenval_mat[k1][k2][0];
		double complex p_2 = S->elastic_eigenval_mat[k1][k2][1];
		double complex p_3 = S->elastic_eigenval_mat[k1][k2][2];
		double complex p_4 = S->elastic_eigenval_mat[k1][k2][3];
		double complex p_5 = S->elastic_eigenval_mat[k1][k2][4];
		double complex p_6 = S->elastic_eigenval_mat[k1][k2][5];
		
		double complex a1_vec_1 = S->elastic_eigenvec_mat[k1][k2][0][0];
		double complex a1_vec_2 = S->elastic_eigenvec_mat[k1][k2][0][1];
		double complex a1_vec_3 = S->elastic_eigenvec_mat[k1][k2][0][2];
		double complex a1_vec_4 = S->elastic_eigenvec_mat[k1][k2][0][3];
		double complex a1_vec_5 = S->elastic_eigenvec_mat[k1][k2][0][4];
		double complex a1_vec_6 = S->elastic_eigenvec_mat[k1][k2][0][5];
			
		double complex a2_vec_1 = S->elastic_eigenvec_mat[k1][k2][1][0];
		double complex a2_vec_2 = S->elastic_eigenvec_mat[k1][k2][1][1];
		double complex a2_vec_3 = S->elastic_eigenvec_mat[k1][k2][1][2];
		double complex a2_vec_4 = S->elastic_eigenvec_mat[k1][k2][1][3];
		double complex a2_vec_5 = S->elastic_eigenvec_mat[k1][k2][1][4];
		double complex a2_vec_6 = S->elastic_eigenvec_mat[k1][k2][1][5];
			
		double complex a3_vec_1 = S->elastic_eigenvec_mat[k1][k2][2][0];
		double complex a3_vec_2 = S->elastic_eigenvec_mat[k1][k2][2][1];
		double complex a3_vec_3 = S->elastic_eigenvec_mat[k1][k2][2][2];
		double complex a3_vec_4 = S->elastic_eigenvec_mat[k1][k2][2][3];
		double complex a3_vec_5 = S->elastic_eigenvec_mat[k1][k2][2][4];
		double complex a3_vec_6 = S->elastic_eigenvec_mat[k1][k2][2][5];
		
		// Calculate z layer
		for(int z = 0; z < C->Nz; z++)
		{
			u_1_B_2Dk_RowMaj[z + C->Nz*(k2 + C->Ny*k1)] = coeffs[0] * a1_vec_1 * cexp(p_1*C->z_axis[z]*I*K_vec) +
														  coeffs[1] * a1_vec_2 * cexp(p_2*C->z_axis[z]*I*K_vec) +
													 	  coeffs[2] * a1_vec_3 * cexp(p_3*C->z_axis[z]*I*K_vec) +
													 	  coeffs[3] * a1_vec_4 * cexp(p_4*C->z_axis[z]*I*K_vec) +
													 	  coeffs[4] * a1_vec_5 * cexp(p_5*C->z_axis[z]*I*K_vec) +
														  coeffs[5] * a1_vec_6 * cexp(p_6*C->z_axis[z]*I*K_vec);

			u_2_B_2Dk_RowMaj[z + C->Nz*(k2 + C->Ny*k1)] = coeffs[0] * a2_vec_1 * cexp(p_1*C->z_axis[z]*I*K_vec) +
														  coeffs[1] * a2_vec_2 * cexp(p_2*C->z_axis[z]*I*K_vec) +
														  coeffs[2] * a2_vec_3 * cexp(p_3*C->z_axis[z]*I*K_vec) +
														  coeffs[3] * a2_vec_4 * cexp(p_4*C->z_axis[z]*I*K_vec) +
														  coeffs[4] * a2_vec_5 * cexp(p_5*C->z_axis[z]*I*K_vec) +
														  coeffs[5] * a2_vec_6 * cexp(p_6*C->z_axis[z]*I*K_vec);

			u_3_B_2Dk_RowMaj[z + C->Nz*(k2 + C->Ny*k1)] = coeffs[0] * a3_vec_1 * cexp(p_1*C->z_axis[z]*I*K_vec) +
														  coeffs[1] * a3_vec_2 * cexp(p_2*C->z_axis[z]*I*K_vec) +
														  coeffs[2] * a3_vec_3 * cexp(p_3*C->z_axis[z]*I*K_vec) +
														  coeffs[3] * a3_vec_4 * cexp(p_4*C->z_axis[z]*I*K_vec) +
														  coeffs[4] * a3_vec_5 * cexp(p_5*C->z_axis[z]*I*K_vec) +
														  coeffs[5] * a3_vec_6 * cexp(p_6*C->z_axis[z]*I*K_vec);


			// du/dx3 
			// analytical expr for du/dx3
			u_1_B_d3_2Dk_RowMaj[z + C->Nz*(k2 + C->Ny*k1)] = (coeffs[0] * a1_vec_1 * cexp(p_1*C->z_axis[z]*I*K_vec)*p_1+
															  coeffs[1] * a1_vec_2 * cexp(p_2*C->z_axis[z]*I*K_vec)*p_2 +
															  coeffs[2] * a1_vec_3 * cexp(p_3*C->z_axis[z]*I*K_vec)*p_3 +
															  coeffs[3] * a1_vec_4 * cexp(p_4*C->z_axis[z]*I*K_vec)*p_4 +
															  coeffs[4] * a1_vec_5 * cexp(p_5*C->z_axis[z]*I*K_vec)*p_5 +
															  coeffs[5] * a1_vec_6 * cexp(p_6*C->z_axis[z]*I*K_vec)*p_6) * I * K_vec;

			u_2_B_d3_2Dk_RowMaj[z + C->Nz*(k2 + C->Ny*k1)] = (coeffs[0] * a2_vec_1 * cexp(p_1*C->z_axis[z]*I*K_vec)*p_1 +
															  coeffs[1] * a2_vec_2 * cexp(p_2*C->z_axis[z]*I*K_vec)*p_2 +
															  coeffs[2] * a2_vec_3 * cexp(p_3*C->z_axis[z]*I*K_vec)*p_3 +
															  coeffs[3] * a2_vec_4 * cexp(p_4*C->z_axis[z]*I*K_vec)*p_4 +
															  coeffs[4] * a2_vec_5 * cexp(p_5*C->z_axis[z]*I*K_vec)*p_5 +
															  coeffs[5] * a2_vec_6 * cexp(p_6*C->z_axis[z]*I*K_vec)*p_6) * I * K_vec;

			u_3_B_d3_2Dk_RowMaj[z + C->Nz*(k2 + C->Ny*k1)] = (coeffs[0] * a3_vec_1 * cexp(p_1*C->z_axis[z]*I*K_vec)*p_1 +
															  coeffs[1] * a3_vec_2 * cexp(p_2*C->z_axis[z]*I*K_vec)*p_2 +
															  coeffs[2] * a3_vec_3 * cexp(p_3*C->z_axis[z]*I*K_vec)*p_3 +
															  coeffs[3] * a3_vec_4 * cexp(p_4*C->z_axis[z]*I*K_vec)*p_4 +
															  coeffs[4] * a3_vec_5 * cexp(p_5*C->z_axis[z]*I*K_vec)*p_5 +
															  coeffs[5] * a3_vec_6 * cexp(p_6*C->z_axis[z]*I*K_vec)*p_6) * I * K_vec;

		 }
		

				
		/* // Construct the solution using PA = LU
		double complex **** elastic_eigenvec_mat; 
		double complex *** elastic_eigenval_mat;
		double complex *** elastic_bc_mat_L_ColMaj;
		double complex *** elastic_bc_mat_U_ColMaj;
		double *** elastic_bc_mat_P_ColMaj;	
		 */
		
	}}

	// kspace origin
	bc[0] = bc_sub_1_2Dk_RowMaj[0];
	bc[1] = bc_sub_2_2Dk_RowMaj[0];
	bc[2] = bc_sub_3_2Dk_RowMaj[0];
	bc[3] = bc_film_1_2Dk_RowMaj[0];
	bc[4] = bc_film_2_2Dk_RowMaj[0];
	bc[5] = bc_film_3_2Dk_RowMaj[0];
	zgemm_(&N, &N, &n6, &n1, &n6, 
		&one, S -> elastic_bc_mat_inv_korigin_ColMaj, &n6, 
		bc, &n6, 
		&zero, coeffs, &n6);
	for(int z = 0; z < C->Nz; z++)
	{
		u_1_B_2Dk_RowMaj[z] = coeffs[0] * C->z_axis[z] + coeffs[1];
		u_2_B_2Dk_RowMaj[z] = coeffs[2] * C->z_axis[z] + coeffs[3];
		u_3_B_2Dk_RowMaj[z] = coeffs[4] * C->z_axis[z] + coeffs[5];
		
		u_1_B_d3_2Dk_RowMaj[z] = coeffs[0];
		u_2_B_d3_2Dk_RowMaj[z] = coeffs[2];
		u_3_B_d3_2Dk_RowMaj[z] = coeffs[4];	
	}		
	
	// IFFT
	// uB is in 2Dk space, each layer of 3D matrix is a 2D fft
	// need to ifft2 slice this
	double complex * u_1_B_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_2_B_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * u_3_B_RowMaj = calloc(C->total, sizeof(double complex));

	IFFT2_slices(C->Nx, C->Ny, C->Nz, u_1_B_2Dk_RowMaj, u_1_B_RowMaj, 0, C->Nz);
	IFFT2_slices(C->Nx, C->Ny, C->Nz, u_2_B_2Dk_RowMaj, u_2_B_RowMaj, 0, C->Nz);
	IFFT2_slices(C->Nx, C->Ny, C->Nz, u_3_B_2Dk_RowMaj, u_3_B_RowMaj, 0, C->Nz);
			
	double complex * e_11_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_12_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_13_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_22_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_23_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_33_B_2Dk_RowMaj = calloc(C->total, sizeof(double complex));

	// Find B 2Dk strain
	for(int i = 0; i < C->total; i++)
	{
		// du1/dx1
		e_11_B_2Dk_RowMaj[i] = I * C->kx_grid_3D_RowMaj[i] * u_1_B_2Dk_RowMaj[i];
		
		// du2/dx2
		e_22_B_2Dk_RowMaj[i] = I * C->ky_grid_3D_RowMaj[i] * u_2_B_2Dk_RowMaj[i];
		
		// du3/dx3 
		e_33_B_2Dk_RowMaj[i] = u_3_B_d3_2Dk_RowMaj[i];
		
		// (du1/dx2 + du2/dx1)*0.5
		e_12_B_2Dk_RowMaj[i] = 0.5 * (I * C->kx_grid_3D_RowMaj[i] * u_2_B_2Dk_RowMaj[i] + 
									  I * C->ky_grid_3D_RowMaj[i] * u_1_B_2Dk_RowMaj[i]);

		// (du1/dx3 + du3/dx1)*0.5
		e_13_B_2Dk_RowMaj[i] = 0.5 * (u_1_B_d3_2Dk_RowMaj[i] + 
									  I * C->kx_grid_3D_RowMaj[i] * u_3_B_2Dk_RowMaj[i]);
		
		// (du2/dx3 + du3/dx2)*0.5
		e_23_B_2Dk_RowMaj[i] = 0.5 * (u_2_B_d3_2Dk_RowMaj[i] + 
									  I * C->ky_grid_3D_RowMaj[i] * u_3_B_2Dk_RowMaj[i]);
	}
		
	double complex * e_11_B_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_12_B_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_13_B_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_22_B_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_23_B_RowMaj = calloc(C->total, sizeof(double complex));
	double complex * e_33_B_RowMaj = calloc(C->total, sizeof(double complex));
	
	IFFT2_slices(C->Nx, C->Ny, C->Nz, e_11_B_2Dk_RowMaj, e_11_B_RowMaj, 0, C->Nz);
	IFFT2_slices(C->Nx, C->Ny, C->Nz, e_12_B_2Dk_RowMaj, e_12_B_RowMaj, 0, C->Nz);
	IFFT2_slices(C->Nx, C->Ny, C->Nz, e_13_B_2Dk_RowMaj, e_13_B_RowMaj, 0, C->Nz);
	IFFT2_slices(C->Nx, C->Ny, C->Nz, e_22_B_2Dk_RowMaj, e_22_B_RowMaj, 0, C->Nz);
	IFFT2_slices(C->Nx, C->Ny, C->Nz, e_23_B_2Dk_RowMaj, e_23_B_RowMaj, 0, C->Nz);
	IFFT2_slices(C->Nx, C->Ny, C->Nz, e_33_B_2Dk_RowMaj, e_33_B_RowMaj, 0, C->Nz);
		
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	// Find A + B
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	for(int i = 0; i < C->total; i++)
	{
		T -> e_11_RowMaj[i] = e_11_A_RowMaj[i] + e_11_B_RowMaj[i]; 
		T -> e_12_RowMaj[i] = e_12_A_RowMaj[i] + e_12_B_RowMaj[i]; 
		T -> e_13_RowMaj[i] = e_13_A_RowMaj[i] + e_13_B_RowMaj[i]; 
		T -> e_22_RowMaj[i] = e_22_A_RowMaj[i] + e_22_B_RowMaj[i]; 
		T -> e_23_RowMaj[i] = e_23_A_RowMaj[i] + e_23_B_RowMaj[i]; 
		T -> e_33_RowMaj[i] = e_33_A_RowMaj[i] + e_33_B_RowMaj[i]; 
		
		T -> u_1_RowMaj[i] = u_1_A_RowMaj[i] + u_1_B_RowMaj[i]; 
		T -> u_2_RowMaj[i] = u_2_A_RowMaj[i] + u_2_B_RowMaj[i]; 
		T -> u_3_RowMaj[i] = u_3_A_RowMaj[i] + u_3_B_RowMaj[i]; 
		T -> TotalStrain_11_RowMaj[i] = e_11_A_RowMaj[i] + e_11_B_RowMaj[i] + C->Us_11; 
		T -> TotalStrain_12_RowMaj[i] = e_12_A_RowMaj[i] + e_12_B_RowMaj[i] + C->Us_12; 
		T -> TotalStrain_13_RowMaj[i] = e_13_A_RowMaj[i] + e_13_B_RowMaj[i]; 
		T -> TotalStrain_22_RowMaj[i] = e_22_A_RowMaj[i] + e_22_B_RowMaj[i] + C->Us_22; 
		T -> TotalStrain_23_RowMaj[i] = e_23_A_RowMaj[i] + e_23_B_RowMaj[i]; 
		T -> TotalStrain_33_RowMaj[i] = e_33_A_RowMaj[i] + e_33_B_RowMaj[i] -
			( C->C12 * C->Us_11 + C->C12 * C->Us_22 ) / C->C44; 

        T -> f1_elastic_RowMaj[i] = 
			-( (C->q11*T->TotalStrain_11_RowMaj[i] + 
			C->q12*T->TotalStrain_22_RowMaj[i] + 
			C->q12*T->TotalStrain_33_RowMaj[i]) * (2*P->P1_RowMaj[i]) + 
			2*C->q44*(T->TotalStrain_12_RowMaj[i] * P->P2_RowMaj[i] + 
			T->TotalStrain_13_RowMaj[i] * P->P3_RowMaj[i]) ) + 
			(4 * (C->b11) * P->P1_RowMaj[i] * P->P1_RowMaj[i] + 
			2 * (C->b12) * ( P->P2_RowMaj[i] * P->P2_RowMaj[i] + 
			P->P3_RowMaj[i] * P->P3_RowMaj[i] )) * P->P1_RowMaj[i];
			
        T -> f2_elastic_RowMaj[i] = 
			-( (C->q11*T->TotalStrain_22_RowMaj[i] + 
			C->q12*T->TotalStrain_11_RowMaj[i] + 
			C->q12*T->TotalStrain_33_RowMaj[i]) * (2*P->P2_RowMaj[i]) + 
			2*C->q44*(T->TotalStrain_12_RowMaj[i] * P->P1_RowMaj[i] + 
			T->TotalStrain_13_RowMaj[i] * P->P3_RowMaj[i]) ) + 
			(4 * (C->b11) * P->P2_RowMaj[i] * P->P2_RowMaj[i] + 
			2 * (C->b12) * ( P->P3_RowMaj[i] * P->P3_RowMaj[i] + 
			P->P1_RowMaj[i] * P->P1_RowMaj[i] )) * P->P2_RowMaj[i];
			
        T -> f3_elastic_RowMaj[i] = 
			-( (C->q11*T->TotalStrain_33_RowMaj[i] +
			C->q12*T->TotalStrain_11_RowMaj[i] + 
			C->q12*T->TotalStrain_22_RowMaj[i]) * (2*P->P3_RowMaj[i]) + 
			2*C->q44*(T->TotalStrain_23_RowMaj[i] * P->P2_RowMaj[i] +
			T->TotalStrain_13_RowMaj[i] * P->P1_RowMaj[i]) ) + 
			(4 * (C->b11) * P->P3_RowMaj[i] * P->P3_RowMaj[i] + 
			2 * (C->b12) * ( P->P1_RowMaj[i] * P->P1_RowMaj[i] + 
			P->P2_RowMaj[i] * P->P2_RowMaj[i] )) * P->P3_RowMaj[i];

	}

	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	// Free memory
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	free(bc_film_1_3Dk_RowMaj);
	free(bc_film_2_3Dk_RowMaj);
	free(bc_film_3_3Dk_RowMaj);
	
	free(bc_film_1_RowMaj);
	free(bc_film_2_RowMaj);
	free(bc_film_3_RowMaj);	
	
	free(bc_film_1_2D_RowMaj);
	free(bc_film_2_2D_RowMaj);
	free(bc_film_3_2D_RowMaj);	
	
	free(bc_sub_1_2D_RowMaj);
	free(bc_sub_2_2D_RowMaj);
	free(bc_sub_3_2D_RowMaj);	

	
	free(bc_film_1_2Dk_RowMaj);
	free(bc_film_2_2Dk_RowMaj);
	free(bc_film_3_2Dk_RowMaj);	
	
	free(bc_sub_1_2Dk_RowMaj);
	free(bc_sub_2_2Dk_RowMaj);
	free(bc_sub_3_2Dk_RowMaj);	
		
	free(u_1_A_3Dk_RowMaj);
	free(u_2_A_3Dk_RowMaj);
	free(u_3_A_3Dk_RowMaj);
	
	free(u_1_A_RowMaj);
	free(u_2_A_RowMaj);
	free(u_3_A_RowMaj);	
	
	free(e_11_A_3Dk_RowMaj);
	free(e_12_A_3Dk_RowMaj);
	free(e_13_A_3Dk_RowMaj);
	free(e_22_A_3Dk_RowMaj);
	free(e_23_A_3Dk_RowMaj);
	free(e_33_A_3Dk_RowMaj);	
	
	free(e_11_A_RowMaj);
	free(e_12_A_RowMaj);
	free(e_13_A_RowMaj);
	free(e_22_A_RowMaj);
	free(e_23_A_RowMaj);
	free(e_33_A_RowMaj);
	
	free(u_1_B_2Dk_RowMaj);
	free(u_2_B_2Dk_RowMaj);
	free(u_3_B_2Dk_RowMaj);
	
	free(u_1_B_d3_2Dk_RowMaj);
	free(u_2_B_d3_2Dk_RowMaj);
	free(u_3_B_d3_2Dk_RowMaj);
	
	free(u_1_B_RowMaj);
	free(u_2_B_RowMaj);
	free(u_3_B_RowMaj);
	
	free(e_11_B_2Dk_RowMaj);
	free(e_12_B_2Dk_RowMaj);
	free(e_13_B_2Dk_RowMaj);
	free(e_22_B_2Dk_RowMaj);
	free(e_23_B_2Dk_RowMaj);
	free(e_33_B_2Dk_RowMaj);
	
	free(e_11_B_RowMaj);
	free(e_12_B_RowMaj);
	free(e_13_B_RowMaj);
	free(e_22_B_RowMaj);
	free(e_23_B_RowMaj);
	free(e_33_B_RowMaj);
}