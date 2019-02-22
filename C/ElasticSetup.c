#include "ElasticSetup.h"

//
// All LAPACK operates on column major 1D arrays!!!!!!
//

// ---------- Elastic Solution Setup ---------- //
void ElasticSetup( struct ConstantsStruct * Constants, 
				   struct ElasticSolutionStruct * ElasticSolution)
{	

	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	// Extract Constants from ConstantsStruct
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	int Nx = Constants -> Nx;
	int Ny = Constants -> Ny;
	int Nz = Constants -> Nz;
	
	double h_film = Constants -> h_film;
	double h_sub = Constants -> h_sub;

	double C[3][3][3][3];
	for(int i = 0; i < 3; i++){
	for(int j = 0; j < 3; j++){
	for(int k = 0; k < 3; k++){
	for(int l = 0; l < 3; l++){
		C[i][j][k][l] = (Constants->C)[i][j][k][l];
	}}}}
	
	double ** kx_grid_2D = Constants -> kx_grid_2D;
	double ** ky_grid_2D = Constants -> ky_grid_2D;
	
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	// Allocate Memory for the Outputs
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	double **** elastic_eigenproblem_mat = AllocLargeArray_4D(Nx, Ny, 6, 6);
	double complex **** elastic_bc_mat = AllocLargeArray_4D_Complex(Nx, Ny, 6, 6);

	double complex **** elastic_bc_mat_inv = AllocLargeArray_4D_Complex(Nx, Ny, 6, 6);
	double complex ** elastic_bc_mat_inv_korigin = AllocLargeArray_2D_Complex(6, 6);
	double complex **** elastic_eigenvec_mat = AllocLargeArray_4D_Complex(Nx, Ny, 6, 6);
	double complex *** elastic_eigenval_mat = AllocLargeArray_3D_Complex(Nx, Ny, 6);
	// double complex *** elastic_scaling_mat = AllocLargeArray_3D_Complex(Nx, Ny, 6);
	double complex **** elastic_bc_mat_L = AllocLargeArray_4D_Complex(Nx, Ny, 6, 6);
	double complex **** elastic_bc_mat_U = AllocLargeArray_4D_Complex(Nx, Ny, 6, 6);
	double **** elastic_bc_mat_P = AllocLargeArray_4D(Nx, Ny, 6, 6);
	
	double complex * Green_11_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	double complex * Green_12_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	double complex * Green_13_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	double complex * Green_21_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	double complex * Green_22_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	double complex * Green_23_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	
	double complex * Green_31_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	double complex * Green_32_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
	double complex * Green_33_RowMaj = calloc(Nx * Ny * Nz, sizeof(double complex));
		
	// -------------------------------------------------------------------------- //
	// -------------------------------------------------------------------------- //
	
	double ** U = AllocLargeArray_2D(3, 3);
	double m[3] = {0, 0, 1};
	for(int i = 0; i < 3; i++){
	for(int k = 0; k < 3; k++){
		for(int j = 0; j < 3; j++){
		for(int l = 0; l < 3; l++){
			U[i][k] = U[i][k] + C[i][j][k][l] * m[j] * m[l];
		}}
	}}
	
	// Find inverse of U
	double * U_inv_1D = calloc(9, sizeof(double));
	Convert2D_2_ColMajor(3, 3, U, U_inv_1D);
	InverseMat(U_inv_1D, 3);

	// Allocate memory for R, W, and mat
	double ** R = AllocLargeArray_2D(3, 3);
	double ** W = AllocLargeArray_2D(3, 3);
	double complex ** mat = AllocLargeArray_2D_Complex(6, 6);
	
	// 1D Column Major arrays
	double * R_1D = calloc(9, sizeof(double));
	double * W_1D = calloc(9, sizeof(double));
	double * N1_1D = calloc(9, sizeof(double));
	double * N2_1D = calloc(9, sizeof(double));
	double * N3_1D = calloc(9, sizeof(double));	
	double * N3_temp_1D = calloc(9, sizeof(double));		
	double complex * mat_1D = calloc(36, sizeof(double complex));

	// N2 = inv(U)'
	for(int i = 0; i < 3; i++){
	for(int j = 0; j < 3; j++){
		N2_1D[j + i * 3] = U_inv_1D[i + j * 3];
	}}
	
	// For LAPACK
	int n3 = 3;
	double neg_one = -1.0;
	double zero = 0.0;
	double one = 1.0;
	char N = 'N';
	char T = 'T';

	// -------------------------------------------------------------------------- //
	// -------------------------------------------------------------------------- //
	//
	// Start Yulan Calculation for Each k1,k2 2D grid point
	//
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	
	double n[3];
	double vec_len;
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
	if(k1 != 0 || k2 != 0) // Origin point
	{
		vec_len = sqrt(kx_grid_2D[k1][k2] * kx_grid_2D[k1][k2] + 
					   ky_grid_2D[k1][k2] * ky_grid_2D[k1][k2]);
		n[0] = kx_grid_2D[k1][k2] / vec_len;
		n[1] = ky_grid_2D[k1][k2] / vec_len;
		n[2] = 0;
		
		// Find R, W
		for(int i = 0; i < 3; i++){
		for(int k = 0; k < 3; k++){
			for(int j = 0; j < 3; j++){
			for(int l = 0; l < 3; l++){
				R[i][k] = R[i][k] + C[i][j][k][l] * n[j] * m[l];
				W[i][k] = W[i][k] + C[i][j][k][l] * n[j] * n[l];
			}}
		}}
		
		// R, W 2D -> R, W 1D Col Maj
		Convert2D_2_ColMajor(3, 3, R, R_1D);
		Convert2D_2_ColMajor(3, 3, W, W_1D);
		
		// N1 = -inv(U) * R'
		dgemm_(&N, &T, &n3, &n3, &n3, 
				&neg_one, U_inv_1D, &n3, 
				R_1D, &n3, 
				&zero, N1_1D, &n3);

		// N3 = R * inv(U) * R' - W
		
		// N3_temp = R * inv(U)
		dgemm_(&N, &N, &n3, &n3, &n3, 
				&one, R_1D, &n3,
				U_inv_1D, &n3, 
				&zero, N3_temp_1D, &n3);
		
		// N3 = N3_temp * R' - W
		memcpy(N3_1D, W_1D, sizeof(double) * 9);
		dgemm_(&N, &T, &n3, &n3, &n3, 
				&one, N3_temp_1D, &n3,
				R_1D, &n3, 
				&neg_one, N3_1D, &n3);
		
		// mat = [N1, N2; N3, N1'];
		for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			mat[i][j] = N1_1D[j * 3 + i];
		}}
		for(int i = 0; i < 3; i++){
		for(int j = 3; j < 6; j++){
			mat[i][j] = N2_1D[(j-3) * 3 + i];
		}}	
		for(int i = 3; i < 6; i++){
		for(int j = 0; j < 3; j++){
			mat[i][j] = N3_1D[j * 3 + (i-3)];
		}}	
		for(int i = 3; i < 6; i++){
		for(int j = 3; j < 6; j++){
			mat[i][j] = N1_1D[(i-3) * 3 + (j-3)];
		}}
		
		for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			elastic_eigenproblem_mat[k1][k2][i][j] = mat[i][j];
		}}
		Convert2D_2_ColMajor_Complex(6, 6, mat, mat_1D);
		
		// Find eigenvalues/eigenvectors of mat
		double complex eigenVectors[36]; double complex eigenValues[6];
		MatrixComplexEigensystem( eigenVectors, eigenValues, mat_1D, 6);
		
		for (int i = 0; i < 6; i++) // Loop through eigenvalues
		{
			elastic_eigenval_mat[k1][k2][i] = eigenValues[i];
			for (int j = 0; j < 6; ++j)
			{
				elastic_eigenvec_mat[k1][k2][j][i] = eigenVectors[j + i * 6];
			}
		}
	
		// Rezero R, W
		for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			R[i][j] = 0;
			W[i][j] = 0;
		}}
		
		// Rezero arrays
		memset(N1_1D, 0, sizeof(double) * 9);
		memset(N3_temp_1D, 0, sizeof(double) * 9);
		memset(N3_1D, 0, sizeof(double) * 9);
	}}}
	
	
	// Free all dynamically allocated (malloc, calloc) memory
	Free_2D_mat(3, 3, U); free(U_inv_1D);
	Free_2D_mat(3, 3, R); Free_2D_mat(3, 3, W);
	free(R_1D); free(W_1D);	free(N1_1D); free(N2_1D);
	free(N3_1D); free(N3_temp_1D);
	Free_2D_mat_Complex(6, 6, mat);	free(mat_1D);
	
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	//
	// Part 2. Construct part of the boundary condition solution matrices
	//
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	
	int n6 = 6; // For LAPACK
	
	// Allocate Memory
	double complex ** bc_mat = AllocLargeArray_2D_Complex(6, 6); // 2D
	
	// 1D Column Major
	double complex * bc_mat_1D = calloc(36, sizeof(double complex));
	double complex * bc_mat_inv_1D = calloc(36, sizeof(double complex));
	double complex * L_1D = calloc(36, sizeof(double complex));
	double complex * U_1D = calloc(36, sizeof(double complex));
	double * P_1D = calloc(36, sizeof(double));

	// Now construct the BC matrices and invert them
	// Invert -> elastic_bc_mat_inv
	// LU decomposition -> elastic_bc_mat_L, U, P_1D
	// A = P * L * U_1D
	//
	// bc_mat * q = boundary conditions (strain, stress BC)
	// q = coefficients for linear combination of the displacement solution
	//
	// elastic_eigenvec_mat[k1][k2][values in eigenvector][eigenvector #]
	
	double complex val;
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
	if(k1 != 0 || k2 != 0)
	{
		vec_len = sqrt(kx_grid_2D[k1][k2] * kx_grid_2D[k1][k2] + 
					   ky_grid_2D[k1][k2] * ky_grid_2D[k1][k2]);
		
		for(int m = 0; m < 6; m++)
		{
			// Scaling, max of eigenvalue * (h_sub : h_film)
			// double complex temp[film_index-sub_index+1];
			// for(int i = 0; i < film_index-sub_index; i++)
			// 	temp[i] = elastic_eigenval_mat[k1][k2][m]*z_axis[sub_index+i];
			// elastic_scaling_mat[k1][k2][m] = MaxValue_Complex(temp, film_index-sub_index+1);
			
			// elastic_scaling_mat[k1][k2][m] = 0;
/* 				val = 
					elastic_eigenvec_mat[k1][k2][k][m] * 
					cexp(elastic_eigenval_mat[k1][k2][m] * h_film * I * vec_len 
							- elastic_scaling_mat[k1][k2][m]) *
					I / vec_len;		 */	
			// Populate the boundary condition matrix
			for(int k = 0; k < 3; k++)
			{
				val = 
					elastic_eigenvec_mat[k1][k2][k][m] * 
					cexp(elastic_eigenval_mat[k1][k2][m] * h_sub * I * vec_len);
				bc_mat[k][m] = val;
			}
			for(int k = 3; k < 6; k++)
			{
				val = 
					elastic_eigenvec_mat[k1][k2][k][m] * 
					cexp(elastic_eigenval_mat[k1][k2][m] * h_film * I * vec_len) *
					I / vec_len;
				bc_mat[k][m] = val;	
			}
		}
		
		for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			elastic_bc_mat[k1][k2][i][j] = bc_mat[i][j];
		}}
		
		// Invert the bc_mat
		// inv(bc_mat)
		Convert2D_2_ColMajor_Complex(6, 6, bc_mat, bc_mat_inv_1D);
		InverseMat_Complex(bc_mat_inv_1D, 6);
		for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			elastic_bc_mat_inv[k1][k2][i][j] = bc_mat_inv_1D[i + j * 6];
		}}

		
		// LU decomposition of bc_mat
		// A = P * L * U 
		Convert2D_2_ColMajor_Complex(6, 6, bc_mat, bc_mat_1D);
		LU_Decomp_Complex(n6, bc_mat_1D, L_1D, U_1D, P_1D);
		for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			elastic_bc_mat_L[k1][k2][i][j] = L_1D[i + j * 6];
			elastic_bc_mat_U[k1][k2][i][j] = U_1D[i + j * 6];
			elastic_bc_mat_P[k1][k2][i][j] = P_1D[i + j * 6];
		}}
		
	}}}
	
	// k1 = k2 = 0 origin point
	double complex korigin_mat[6][6] = { {h_sub,          1,             0,  0,             0,  0},
										 {0,              0,         h_sub,  1,             0,  0},
										 {0,              0,             0,  0,         h_sub,  1},
										 {C[0][2][0][2],  0, C[0][2][1][2],  0, C[0][2][2][2],  0},
										 {C[1][2][0][2],  0, C[1][2][1][2],  0, C[1][2][2][2],  0},
										 {C[2][2][0][2],  0, C[2][2][1][2],  0, C[2][2][2][2],  0} };
	double complex * korigin_mat_1D = calloc(36, sizeof(double complex));
	for(int i = 0; i < 6; i++){
	for(int j = 0; j < 6; j++){
		korigin_mat_1D[i + j * 6] = korigin_mat[i][j];
	}}
	
	InverseMat_Complex(korigin_mat_1D, 6);
	for(int i = 0; i < 6; i++){
	for(int j = 0; j < 6; j++){
		elastic_bc_mat_inv_korigin[i][j] = 
			korigin_mat_1D[i + j * 6];
	}}

	// Free all dynamically allocated memory
	free(bc_mat_inv_1D); free(korigin_mat_1D);
	Free_2D_mat_Complex(6, 6, bc_mat);
	free(bc_mat_1D);
	free(L_1D); free(U_1D); free(P_1D);
		
	// ----- Find Green's Tensor ----- //
	double K_vec[3] = {0, 0, 0};
	int index_RowMaj = 0;
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
	for(int k3 = 0; k3 < Nz; k3++){
	if(k1 != 0 || k2 != 0 || k3 != 0){
		
		index_RowMaj = k3 + Nz * (k2 + Ny * k1);
		
		K_vec[0] = Constants -> kx_grid_3D_RowMaj[index_RowMaj];
		K_vec[1] = Constants -> ky_grid_3D_RowMaj[index_RowMaj];
		K_vec[2] = Constants -> kz_grid_3D_RowMaj[index_RowMaj];
		
		// column major
		double * g_inv = calloc(9, sizeof(double));
		double * g = calloc(9, sizeof(double));
		int index_ColMaj = 0;
		
		for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			for(int l = 0; l < 3; l++){
			for(int k = 0; k < 3; k++){
				
				index_ColMaj = i + 3 * j;
				g_inv[index_ColMaj] = g_inv[index_ColMaj] + 
					Constants->C[i][k][j][l] * K_vec[k] * K_vec[l];
			
			}}
		}}
		
		// Find inv(g_inv)
		memcpy(g, g_inv, sizeof(double) * 9);
		InverseMat(g, 3);

		// g is col maj
		Green_11_RowMaj[index_RowMaj] = g[0];
		Green_12_RowMaj[index_RowMaj] = g[3];
		Green_13_RowMaj[index_RowMaj] = g[6];

		Green_21_RowMaj[index_RowMaj] = g[1];
		Green_22_RowMaj[index_RowMaj] = g[4];
		Green_23_RowMaj[index_RowMaj] = g[7];

		Green_31_RowMaj[index_RowMaj] = g[2];
		Green_32_RowMaj[index_RowMaj] = g[5];
		Green_33_RowMaj[index_RowMaj] = g[8];
		
		free(g); free(g_inv);
				
	}}}}
	
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//

	// Convert to column major
	// 		elastic_bc_mat_inv
	// 		elastic_bc_mat_inv_korigin
	// 		elastic_bc_mat_L
	// 		elastic_bc_mat_U
	// 		elastic_bc_mat_P
	double complex *** elastic_bc_mat_inv_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 36);
	double complex * elastic_bc_mat_inv_korigin_ColMaj = calloc(36, sizeof(double complex));
	double complex *** elastic_bc_mat_L_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 36);
	double complex *** elastic_bc_mat_L_inv_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 36);
	double complex *** elastic_bc_mat_U_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 36);
	double complex *** elastic_bc_mat_P_ColMaj = AllocLargeArray_3D_Complex(Nx, Ny, 36);

	Convert2D_2_ColMajor_Complex(6, 6, 
		elastic_bc_mat_inv_korigin, 
		elastic_bc_mat_inv_korigin_ColMaj);
			
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
		Convert2D_2_ColMajor_Complex(6, 6, 
			elastic_bc_mat_inv[k1][k2], 
			elastic_bc_mat_inv_ColMaj[k1][k2]);
			
		Convert2D_2_ColMajor_Complex(6, 6, 
			elastic_bc_mat_L[k1][k2], 
			elastic_bc_mat_L_ColMaj[k1][k2]);
		Convert2D_2_ColMajor_Complex(6, 6, 
			elastic_bc_mat_U[k1][k2], 
			elastic_bc_mat_U_ColMaj[k1][k2]);
		Convert2D_2_ColMajor_Double_2_Complex(6, 6, 
			elastic_bc_mat_P[k1][k2], 
			elastic_bc_mat_P_ColMaj[k1][k2]);		
		
	}}
	
	// Find elastic_bc_mat_L_inv_ColMaj
	double complex temp[36];
	for(int k1 = 0; k1 < Nx; k1++){
	for(int k2 = 0; k2 < Ny; k2++){
		memcpy(temp, elastic_bc_mat_L_ColMaj[k1][k2], sizeof(double complex) * 36);
		InverseMat_Complex(temp, 6);
		memcpy(elastic_bc_mat_L_inv_ColMaj[k1][k2], temp, sizeof(double complex) * 36);
	}}	
	
	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//	
	Free_4D_mat_Complex(Nx, Ny, 6, 6, elastic_bc_mat_inv);
	Free_2D_mat_Complex(6, 6, elastic_bc_mat_inv_korigin);
	Free_4D_mat_Complex(Nx, Ny, 6, 6, elastic_bc_mat_L);
	Free_4D_mat_Complex(Nx, Ny, 6, 6, elastic_bc_mat_U);
	Free_4D_mat(Nx, Ny, 6, 6, elastic_bc_mat_P);

	// --------------------------------------------------------------------------//
	// --------------------------------------------------------------------------//
	// Outputs
	ElasticSolution -> elastic_eigenproblem_mat = elastic_eigenproblem_mat;
	ElasticSolution -> elastic_bc_mat = elastic_bc_mat;
	
	//ElasticSolution -> elastic_scaling_mat = elastic_scaling_mat;
	ElasticSolution -> elastic_eigenval_mat = elastic_eigenval_mat;
	ElasticSolution -> elastic_eigenvec_mat = elastic_eigenvec_mat;
	
	ElasticSolution -> elastic_bc_mat_inv_ColMaj = elastic_bc_mat_inv_ColMaj;
	ElasticSolution -> elastic_bc_mat_inv_korigin_ColMaj = elastic_bc_mat_inv_korigin_ColMaj;
	ElasticSolution -> elastic_bc_mat_L_ColMaj = elastic_bc_mat_L_ColMaj;
	ElasticSolution -> elastic_bc_mat_L_inv_ColMaj = elastic_bc_mat_L_inv_ColMaj;
	ElasticSolution -> elastic_bc_mat_U_ColMaj = elastic_bc_mat_U_ColMaj;
	ElasticSolution -> elastic_bc_mat_P_ColMaj = elastic_bc_mat_P_ColMaj;
	
	ElasticSolution -> Green_11_RowMaj = Green_11_RowMaj;
	ElasticSolution -> Green_12_RowMaj = Green_12_RowMaj;
	ElasticSolution -> Green_13_RowMaj = Green_13_RowMaj;

	ElasticSolution -> Green_21_RowMaj = Green_21_RowMaj;
	ElasticSolution -> Green_22_RowMaj = Green_22_RowMaj;
	ElasticSolution -> Green_23_RowMaj = Green_23_RowMaj;

	ElasticSolution -> Green_31_RowMaj = Green_31_RowMaj;
	ElasticSolution -> Green_32_RowMaj = Green_32_RowMaj;
	ElasticSolution -> Green_33_RowMaj = Green_33_RowMaj;
	
}


void FreeElasticSolution( int Nx, int Ny, struct ElasticSolutionStruct * ElasticSolution)
{
	Free_3D_mat_Complex(Nx, Ny, 6, ElasticSolution -> elastic_eigenval_mat);
	Free_4D_mat_Complex(Nx, Ny, 6, 6, ElasticSolution -> elastic_eigenvec_mat);
	
	Free_4D_mat(Nx, Ny, 6, 6, ElasticSolution -> elastic_eigenproblem_mat);
	Free_4D_mat_Complex(Nx, Ny, 6, 6, ElasticSolution -> elastic_bc_mat);	
	
	//(Nx, Ny, 6, ElasticSolution -> elastic_scaling_mat);
    Free_3D_mat_Complex(Nx, Ny, 36, ElasticSolution -> elastic_bc_mat_inv_ColMaj);
    free(ElasticSolution -> elastic_bc_mat_inv_korigin_ColMaj);
    Free_3D_mat_Complex(Nx, Ny, 36, ElasticSolution -> elastic_bc_mat_L_ColMaj);
    Free_3D_mat_Complex(Nx, Ny, 36, ElasticSolution -> elastic_bc_mat_L_inv_ColMaj);
    Free_3D_mat_Complex(Nx, Ny, 36, ElasticSolution -> elastic_bc_mat_U_ColMaj);
    Free_3D_mat_Complex(Nx, Ny, 36, ElasticSolution -> elastic_bc_mat_P_ColMaj);
	
	free(ElasticSolution->Green_11_RowMaj);
	free(ElasticSolution->Green_12_RowMaj);
	free(ElasticSolution->Green_13_RowMaj);

	free(ElasticSolution->Green_21_RowMaj);
	free(ElasticSolution->Green_22_RowMaj);
	free(ElasticSolution->Green_23_RowMaj);

	free(ElasticSolution->Green_31_RowMaj);
	free(ElasticSolution->Green_32_RowMaj);
	free(ElasticSolution->Green_33_RowMaj);
	
}