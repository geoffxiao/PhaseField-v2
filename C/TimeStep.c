#include "TimeStep.h"

// Find time step 1 from time step 0
void TimeStep_1stOrd(struct PVector * P_prev, 
	struct EnergiesStruct * F_prev,
	struct ConstantsStruct * C, 
	struct TimeSolutionStruct * S, 
	struct PVector * P)
{
	// Load Constants
	double dt = C -> dt;
	int Nx = C -> Nx;
	int Ny = C -> Ny;
	int Nz = C -> Nz;
	int total = C -> total;
	double G11 = C -> G11;
	double H44 = C -> H44;
	
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Particular Solution
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	double complex * P1_A_3Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_A_3Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_A_3Dk_RowMaj = calloc(total, sizeof(double complex));

	for(int i = 0; i < total; i++)
	{
		P1_A_3Dk_RowMaj[i] = (( P_prev->P1_3Dk_RowMaj[i] + dt * -F_prev->f1_3Dk_RowMaj[i] ) / 
					( 1.0 + dt * ( H44 * (C->ky_grid_3D_RowMaj[i]*C->ky_grid_3D_RowMaj[i] 
											  + C->kz_grid_3D_RowMaj[i]*C->kz_grid_3D_RowMaj[i] ) + 
								(G11 * C->kx_grid_3D_RowMaj[i] * C->kx_grid_3D_RowMaj[i]) ) ));
	 
		P2_A_3Dk_RowMaj[i] = (( P_prev->P2_3Dk_RowMaj[i] + dt * -F_prev->f2_3Dk_RowMaj[i] ) / 
					( 1.0 + dt * ( H44 * (C->kx_grid_3D_RowMaj[i]*C->kx_grid_3D_RowMaj[i] 
											  + C->kz_grid_3D_RowMaj[i]*C->kz_grid_3D_RowMaj[i] ) 
								+ (G11 * C->ky_grid_3D_RowMaj[i] * C->ky_grid_3D_RowMaj[i]) ) ));
	 
		P3_A_3Dk_RowMaj[i] = (( P_prev->P3_3Dk_RowMaj[i] + dt * -F_prev->f3_3Dk_RowMaj[i] ) / 
					( 1.0 + dt * ( H44 * (C->kx_grid_3D_RowMaj[i]*C->kx_grid_3D_RowMaj[i] 
											  + C->ky_grid_3D_RowMaj[i]*C->ky_grid_3D_RowMaj[i] ) 
								+ (G11 * C->kz_grid_3D_RowMaj[i] * C->kz_grid_3D_RowMaj[i]) ) ));
	}
	
	double complex * P1_A_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_A_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_A_RowMaj = calloc(total, sizeof(double complex));
	
	IFFT3(Nx, Ny, Nz, P1_A_3Dk_RowMaj, P1_A_RowMaj);
	IFFT3(Nx, Ny, Nz, P2_A_3Dk_RowMaj, P2_A_RowMaj);
	IFFT3(Nx, Ny, Nz, P3_A_3Dk_RowMaj, P3_A_RowMaj);
		
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Get film and interface boundary conditions
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	double complex * bc_film_1_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_2_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_3_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));

	double complex * bc_int_1_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_2_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_3_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
				
	// extract the z = film_index and z = sub_index points in the 2D plane
	int index_RowMajor_3D = 0;
	int index_RowMajor_2D = 0;
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		
		index_RowMajor_3D = C->film_index + Nz * (j + Ny * i);
		index_RowMajor_2D = j + Ny * i;
		
		bc_film_1_2D_RowMaj[index_RowMajor_2D] = -P1_A_RowMaj[index_RowMajor_3D];
		bc_film_2_2D_RowMaj[index_RowMajor_2D] = -P2_A_RowMaj[index_RowMajor_3D];
		bc_film_3_2D_RowMaj[index_RowMajor_2D] = -P3_A_RowMaj[index_RowMajor_3D];
		
		index_RowMajor_3D = C->interface_index + Nz * (j + Ny * i);		
		bc_int_1_2D_RowMaj[index_RowMajor_2D] = -P1_A_RowMaj[index_RowMajor_3D];
		bc_int_2_2D_RowMaj[index_RowMajor_2D] = -P2_A_RowMaj[index_RowMajor_3D];
		bc_int_3_2D_RowMaj[index_RowMajor_2D] = -P3_A_RowMaj[index_RowMajor_3D];
		
	}}
	
	
	// FFT2 the bc
	double complex * bc_film_1_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_2_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_3_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));

	double complex * bc_int_1_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_2_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_3_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
				
	FFT2(Nx, Ny, bc_film_1_2D_RowMaj, bc_film_1_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_film_2_2D_RowMaj, bc_film_2_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_film_3_2D_RowMaj, bc_film_3_2Dk_RowMaj); 

	FFT2(Nx, Ny, bc_int_1_2D_RowMaj, bc_int_1_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_int_2_2D_RowMaj, bc_int_2_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_int_3_2D_RowMaj, bc_int_3_2Dk_RowMaj); 

	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Homo solution
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //	
	double complex * P1_B_2Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_B_2Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_B_2Dk_RowMaj = calloc(total, sizeof(double complex));
	
	double complex bc1[2] = {0, 0}; 
	double complex bc2[2] = {0, 0}; 
	double complex bc3[2] = {0, 0}; 
	int index_RowMaj; 
	double complex coeffs1[2] = {0, 0};
	double complex coeffs2[2] = {0, 0};
	double complex coeffs3[2] = {0, 0};
	double complex one = 1.0; double complex zero = 0.0;
	char N = 'N';
	int n2 = 2;
	int n1 = 1;	
	
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){

		index_RowMaj = k2 + C -> Ny * k1;
		
		bc1[0] = bc_int_1_2Dk_RowMaj[index_RowMaj];
		bc1[1] = bc_film_1_2Dk_RowMaj[index_RowMaj];

		zgemm_(&N, &N, &n2, &n1, &n2, 
				&one, S -> time_bc_mat_inv_1_1storder_ColMaj[k1][k2], &n2, 
				bc1, &n2, 
				&zero, coeffs1, &n2);

		bc2[0] = bc_int_2_2Dk_RowMaj[index_RowMaj];
		bc2[1] = bc_film_2_2Dk_RowMaj[index_RowMaj];

		zgemm_(&N, &N, &n2, &n1, &n2, 
				&one, S -> time_bc_mat_inv_2_1storder_ColMaj[k1][k2], &n2, 
				bc2, &n2, 
				&zero, coeffs2, &n2);
				
		bc3[0] = bc_int_3_2Dk_RowMaj[index_RowMaj];
		bc3[1] = bc_film_3_2Dk_RowMaj[index_RowMaj];

		zgemm_(&N, &N, &n2, &n1, &n2, 
				&one, S -> time_bc_mat_inv_3_1storder_ColMaj[k1][k2], &n2, 
				bc3, &n2, 
				&zero, coeffs3, &n2);	

		// Fill in the z direction
		for(int z = 0; z < Nz; z++)
		{
			P1_B_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				coeffs1[0] * cexp(S->time_q_1_1storder[k1][k2] * C->z_axis[z]) + 
				coeffs1[1] * cexp(-S->time_q_1_1storder[k1][k2] * C->z_axis[z]);
			P2_B_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				coeffs2[0] * cexp(S->time_q_2_1storder[k1][k2] * C->z_axis[z]) + 
				coeffs2[1] * cexp(-S->time_q_2_1storder[k1][k2] * C->z_axis[z]);
			P3_B_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				coeffs3[0] * cexp(S->time_q_3_1storder[k1][k2] * C->z_axis[z]) + 
				coeffs3[1] * cexp(-S->time_q_3_1storder[k1][k2] * C->z_axis[z]);

		}			
	}}
	
	double complex * P1_B_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_B_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_B_RowMaj = calloc(total, sizeof(double complex));
	
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, P1_B_2Dk_RowMaj, P1_B_RowMaj, 0, C -> Nz);
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, P2_B_2Dk_RowMaj, P2_B_RowMaj, 0, C -> Nz);
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, P3_B_2Dk_RowMaj, P3_B_RowMaj, 0, C -> Nz);	
	
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Total
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	for(int i = 0; i < total; i++)
	{
		P -> P1_RowMaj[i] = (P1_B_RowMaj[i] + P1_A_RowMaj[i]) * C->in_film_RowMaj[i];
		P -> P2_RowMaj[i] = (P2_B_RowMaj[i] + P2_A_RowMaj[i]) * C->in_film_RowMaj[i];
		P -> P3_RowMaj[i] = (P3_B_RowMaj[i] + P3_A_RowMaj[i]) * C->in_film_RowMaj[i];
	}
	
	// FFT3
	FFT3(Nx, Ny, Nz, P -> P1_RowMaj, P -> P1_3Dk_RowMaj);
	FFT3(Nx, Ny, Nz, P -> P2_RowMaj, P -> P2_3Dk_RowMaj);
	FFT3(Nx, Ny, Nz, P -> P3_RowMaj, P -> P3_3Dk_RowMaj);

	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Free
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //	
	free(bc_film_1_2D_RowMaj);
	free(bc_film_2_2D_RowMaj);
	free(bc_film_3_2D_RowMaj);
	
	free(bc_int_1_2D_RowMaj);
	free(bc_int_2_2D_RowMaj);
	free(bc_int_3_2D_RowMaj);
	
	free(P1_A_3Dk_RowMaj);
	free(P2_A_3Dk_RowMaj);
	free(P3_A_3Dk_RowMaj);
	
	free(P1_A_RowMaj);
	free(P2_A_RowMaj);
	free(P3_A_RowMaj);	
	
	free(P1_B_2Dk_RowMaj);
	free(P2_B_2Dk_RowMaj);
	free(P3_B_2Dk_RowMaj);
	
	free(P1_B_RowMaj);
	free(P2_B_RowMaj);
	free(P3_B_RowMaj);	
	
	free(bc_film_1_2Dk_RowMaj);
	free(bc_film_2_2Dk_RowMaj);
	free(bc_film_3_2Dk_RowMaj);
	
	free(bc_int_1_2Dk_RowMaj);
	free(bc_int_2_2Dk_RowMaj);
	free(bc_int_3_2Dk_RowMaj);	
}








/* // Find time step 2 from time step 0 and step 1
void TimeStep_2ndOrd(struct PVector * P_0, struct PVector * P_1,
	struct EnergiesStruct * F_0, struct EnergiesStruct * F_1,
	struct ConstantsStruct * C, 
	struct TimeSolutionStruct * S, 
	struct PVector * P)
{
	// Load Constants
	double dt = C -> dt;
	int Nx = C -> Nx;
	int Ny = C -> Ny;
	int Nz = C -> Nz;
	int interface_index = C -> interface_index;
	double h_int = C -> h_int;
	int film_index = C -> film_index;
	double h_film = C -> h_film;
	double G11 = C -> G11;
	double H44 = C -> H44;
	
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Particular Solution
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	double complex * P1_A_3Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_A_3Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_A_3Dk_RowMaj = calloc(total, sizeof(double complex));

	for(int i = 0; i < total; i++)
	{
		P1_A_3Dk_RowMaj[i] = (( P_prev->P1_3Dk_RowMaj[i] + C->dt * -F_prev->f1_3Dk_RowMaj[i] ) / 
					( 1.0 + dt * ( H44 * (C->ky_grid_3D_RowMaj[i]*C->ky_grid_3D_RowMaj[i] 
											  + C->kz_grid_3D_RowMaj[i]*C->kz_grid_3D_RowMaj[i] ) + 
								(G11 * C->kx_grid_3D_RowMaj[i] * C->kx_grid_3D_RowMaj[i]) ) ));
	 
		P2_A_3Dk_RowMaj[i] = (( P_prev->P2_3Dk_RowMaj[i] + dt * -F_prev->f2_3Dk_RowMaj[i] ) / 
					( 1.0 + dt * ( H44 * (C->kx_grid_3D_RowMaj[i]*C->kx_grid_3D_RowMaj[i] 
											  + C->kz_grid_3D_RowMaj[i]*C->kz_grid_3D_RowMaj[i] ) 
								+ (G11 * C->ky_grid_3D_RowMaj[i] * C->ky_grid_3D_RowMaj[i]) ) ));
	 
		P3_A_3Dk_RowMaj[i] = (( P_prev->P3_3Dk_RowMaj[i] + dt * -F_prev->f3_3Dk_RowMaj[i] ) / 
					( 1.0 + dt * ( H44 * (C->kx_grid_3D_RowMaj[i]*C->kx_grid_3D_RowMaj[i] 
											  + C->ky_grid_3D_RowMaj[i]*C->ky_grid_3D_RowMaj[i] ) 
								+ (G11 * C->ky_grid_3D_RowMaj[i] * C->ky_grid_3D_RowMaj[i]) ) ));
	}
	
	double complex * P1_A_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_A_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_A_RowMaj = calloc(total, sizeof(double complex));
	
	IFFT3(Nx, Ny, Nz, P1_A_3Dk_RowMaj, P1_A_RowMaj);
	IFFT3(Nx, Ny, Nz, P2_A_3Dk_RowMaj, P2_A_RowMaj);
	IFFT3(Nx, Ny, Nz, P3_A_3Dk_RowMaj, P3_A_RowMaj);
		
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Get film and interface boundary conditions
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	double complex * bc_film_1_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_2_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_3_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));

	double complex * bc_int_1_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_2_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_3_2D_RowMaj = calloc(Nx * Ny, sizeof(double complex));
				
	// extract the z = film_index and z = sub_index points in the 2D plane
	int index_RowMajor_3D = 0;
	int index_RowMajor_2D = 0;
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
		
		index_RowMajor_3D = C->film_index + Nz * (j + Ny * i);
		index_RowMajor_2D = j + Ny * i;
		
		bc_film_1_2D_RowMaj[index_RowMajor_2D] = -P1_A_RowMaj[index_RowMajor_3D];
		bc_film_2_2D_RowMaj[index_RowMajor_2D] = -P2_A_RowMaj[index_RowMajor_3D];
		bc_film_3_2D_RowMaj[index_RowMajor_2D] = -P3_A_RowMaj[index_RowMajor_3D];
		
		index_RowMajor_3D = C->interface_index + Nz * (j + Ny * i);		
		bc_int_1_2D_RowMaj[index_RowMajor_2D] = -P1_A_RowMaj[index_RowMajor_3D];
		bc_int_2_2D_RowMaj[index_RowMajor_2D] = -P2_A_RowMaj[index_RowMajor_3D];
		bc_int_3_2D_RowMaj[index_RowMajor_2D] = -P3_A_RowMaj[index_RowMajor_3D];
		
	}}
	
	
	// FFT2 the bc
	double complex * bc_film_1_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_2_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_film_3_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));

	double complex * bc_int_1_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_2_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
	double complex * bc_int_3_2Dk_RowMaj = calloc(Nx * Ny, sizeof(double complex));
				
	FFT2(Nx, Ny, bc_film_1_2D_RowMaj, bc_film_1_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_film_2_2D_RowMaj, bc_film_2_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_film_3_2D_RowMaj, bc_film_3_2Dk_RowMaj); 

	FFT2(Nx, Ny, bc_int_1_2D_RowMaj, bc_int_1_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_int_2_2D_RowMaj, bc_int_2_2Dk_RowMaj); 
	FFT2(Nx, Ny, bc_int_3_2D_RowMaj, bc_int_3_2Dk_RowMaj); 

	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Homo solution
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //	
	double complex * P1_B_2Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_B_2Dk_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_B_2Dk_RowMaj = calloc(total, sizeof(double complex));
	
	double complex bc[2] = {0, 0}; 
	int index_RowMaj; 
	double complex coeffs1[2] = {0, 0};
	double complex coeffs2[2] = {0, 0};
	double complex coeffs3[2] = {0, 0};
	double complex one = 1.0; double complex zero = 0.0;
	char N = 'N';
	int n2 = 2;
	int n1 = 1;	
	
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){

		index_RowMaj = k2 + C -> Ny * k1;
		
		bc[0] = bc_int_1_2Dk_RowMaj[index_RowMaj];
		bc[1] = bc_film_1_2Dk_RowMaj[index_RowMaj];

		zgemm_(&N, &N, &n2, &n1, &n2, 
				&one, S -> time_bc_mat_inv_1_1storder_ColMaj[k1][k2], &n2, 
				bc, &n2, 
				&zero, coeffs1, &n2);

		bc[0] = bc_int_2_2Dk_RowMaj[index_RowMaj];
		bc[1] = bc_film_2_2Dk_RowMaj[index_RowMaj];

		zgemm_(&N, &N, &n2, &n1, &n2, 
				&one, S -> time_bc_mat_inv_2_1storder_ColMaj[k1][k2], &n2, 
				bc, &n2, 
				&zero, coeffs2, &n2);
				
		bc[0] = bc_int_3_2Dk_RowMaj[index_RowMaj];
		bc[1] = bc_film_3_2Dk_RowMaj[index_RowMaj];

		zgemm_(&N, &N, &n2, &n1, &n2, 
				&one, S -> time_bc_mat_inv_3_1storder_ColMaj[k1][k2], &n2, 
				bc, &n2, 
				&zero, coeffs3, &n2);	

		// Fill in the z direction
		for(int z = 0; z < Nz; z++)
		{
			P1_B_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				coeffs1[0] * cexp(S->time_q_1_1storder[k1][k2] * C->z_axis[z]) + 
				coeffs1[1] * cexp(-S->time_q_1_1storder[k1][k2] * C->z_axis[z]);
			P2_B_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				coeffs2[0] * cexp(S->time_q_2_1storder[k1][k2] * C->z_axis[z]) + 
				coeffs2[1] * cexp(-S->time_q_2_1storder[k1][k2] * C->z_axis[z]);
			P3_B_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				coeffs3[0] * cexp(S->time_q_3_1storder[k1][k2] * C->z_axis[z]) + 
				coeffs3[1] * cexp(-S->time_q_3_1storder[k1][k2] * C->z_axis[z]);

		}			
	}}
	
	double complex * P1_B_RowMaj = calloc(total, sizeof(double complex));
	double complex * P2_B_RowMaj = calloc(total, sizeof(double complex));
	double complex * P3_B_RowMaj = calloc(total, sizeof(double complex));
	
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, P1_B_2Dk_RowMaj, P1_B_RowMaj, 0, C -> Nz);
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, P2_B_2Dk_RowMaj, P2_B_RowMaj, 0, C -> Nz);
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, P3_B_2Dk_RowMaj, P3_B_RowMaj, 0, C -> Nz);	
	
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	// Total
	// ------------------------------------------------------------------- //
	// ------------------------------------------------------------------- //
	for(int i = 0; i < total; i++)
	{
		P -> P1_RowMaj[i] = (P1_B_RowMaj[i] + P1_A_RowMaj[i]) ; C->in_film_RowMaj[i];
		P -> P2_RowMaj[i] = (P2_B_RowMaj[i] + P2_A_RowMaj[i]) ; C->in_film_RowMaj[i];
		P -> P3_RowMaj[i] = (P3_B_RowMaj[i] + P3_A_RowMaj[i]) ; C->in_film_RowMaj[i];
	}
	
	// FFT3
	FFT3(Nx, Ny, Nz, P -> P1_RowMaj, P -> P1_3Dk_RowMaj);
	FFT3(Nx, Ny, Nz, P -> P2_RowMaj, P -> P2_3Dk_RowMaj);
	FFT3(Nx, Ny, Nz, P -> P3_RowMaj, P -> P3_3Dk_RowMaj);
	
	for(int i = 0; i < Nx*Ny*Nz; i++)
		printf("%g + %g * 1i\n", creal(P -> P3_RowMaj[i]), 
				cimag(P -> P3_RowMaj[i]));
	
	// Free
	free(bc_film_1_2D_RowMaj);
	free(bc_film_2_2D_RowMaj);
	free(bc_film_3_2D_RowMaj);
	
	free(bc_int_1_2D_RowMaj);
	free(bc_int_2_2D_RowMaj);
	free(bc_int_3_2D_RowMaj);
	
	free(P1_A_3Dk_RowMaj);
	free(P2_A_3Dk_RowMaj);
	free(P3_A_3Dk_RowMaj);
	
	free(P1_A_RowMaj);
	free(P2_A_RowMaj);
	free(P3_A_RowMaj);	
	
	free(P1_B_2Dk_RowMaj);
	free(P2_B_2Dk_RowMaj);
	free(P3_B_2Dk_RowMaj);
	
	free(P1_B_RowMaj);
	free(P2_B_RowMaj);
	free(P3_B_RowMaj);	
}
 */