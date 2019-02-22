#include "CalcElecEnergy.h"

// Find particular solution
void ParticularSolution(struct ConstantsStruct * C, struct PVector * P, 
	struct ElectricWorkspaceStruct * E,
	int start, int end)
{
	for(int i = start; i < end; i++)
	{
		E -> potential_A_3Dk_RowMaj[i] = 
			-1 * I * (C -> kx_grid_3D_RowMaj[i] * P -> P1_3Dk_RowMaj[i] + 
					  C -> ky_grid_3D_RowMaj[i] * P -> P2_3Dk_RowMaj[i] + 
					  C -> kz_grid_3D_RowMaj[i] * P -> P3_3Dk_RowMaj[i]) /
            ( C -> permittivity_0 * 
				(C -> k_electric_11 * C -> kx_grid_3D_RowMaj[i] * C -> kx_grid_3D_RowMaj[i] + 
				 C -> k_electric_22 * C -> ky_grid_3D_RowMaj[i] * C -> ky_grid_3D_RowMaj[i] + 
				 C -> k_electric_33 * C -> kz_grid_3D_RowMaj[i] * C -> kz_grid_3D_RowMaj[i]) );

		E -> E_A_1_3Dk_RowMaj[i] = -1 * I * C -> kx_grid_3D_RowMaj[i] * E -> potential_A_3Dk_RowMaj[i];
		E -> E_A_2_3Dk_RowMaj[i] = -1 * I * C -> ky_grid_3D_RowMaj[i] * E -> potential_A_3Dk_RowMaj[i];
		E -> E_A_3_3Dk_RowMaj[i] = -1 * I * C -> kz_grid_3D_RowMaj[i] * E -> potential_A_3Dk_RowMaj[i];
	}
}
// ---------------------------------------------------------------------------//
// ---------------------------------------------------------------------------//	


void HomoSolution(struct ConstantsStruct * C, struct ElectricWorkspaceStruct * E, 
	struct ElectricSolutionStruct * S, 
	int start, int end)
{	
	double complex bc[2] = {0, 0}; 
	int index_RowMaj; 
	double complex coeffs[2] = {0, 0};
	double complex one = 1.0; double complex zero = 0.0;
	char N = 'N';
	int n2 = 2;
	int n1 = 1;
	// Calculate coefficients
	for(int k1 = start; k1 < end; k1++){
	for(int k2 = 0; k2 < C -> Ny; k2++){
	if(k1 != 0 || k2 != 0){
		
		// k1, k2 -> Row Major index
		index_RowMaj = k2 + C -> Ny * k1;
		
		// RHS [interface; film]
		bc[0] = E -> potential_bc_interface_2Dk_RowMaj[index_RowMaj];
		bc[1] = E -> potential_bc_film_2Dk_RowMaj[index_RowMaj];

		// bc_mat_inv * RHS 
		// A = m x k
		// B = k x n
		// C = m x n
		// Inputs: m, n, k
		//	       m = 2, n = 1, k = 2
		zgemm_(&N, &N, &n2, &n1, &n2, 
				&one, S -> electric_bc_mat_inv_ColMaj[k1][k2], &n2, 
				bc, &n2, 
				&zero, coeffs, &n2);
				
		for(int z = 0; z < C -> Nz; z++)
		{
			E -> potential_B_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				coeffs[0] * cexp(  C -> z_axis[z]  * S -> electric_p_mat[k1][k2]) + 
				coeffs[1] * cexp(-(C -> z_axis[z]) * S -> electric_p_mat[k1][k2]);
			E -> potential_B_d3_2Dk_RowMaj[index_RowMaj * C -> Nz + z] = 
				S -> electric_p_mat[k1][k2] * coeffs[0] * cexp(  C -> z_axis[z]  * S -> electric_p_mat[k1][k2]) -
				S -> electric_p_mat[k1][k2] * coeffs[1] * cexp(-(C -> z_axis[z]) * S -> electric_p_mat[k1][k2]);
		}
	}}}	
}
// ---------------------------------------------------------------------------//
// ---------------------------------------------------------------------------//	







//
// Mulithread???
//

// ---------------------------------------------------------------------------//
// ---------------------------------------------------------------------------//

void CalcElecEnergy(struct PVector * P, struct ConstantsStruct * C, 
	struct ElectricSolutionStruct * S, struct ElectricWorkspaceStruct * E,
	struct TotalStateStruct * T)
{
/* 	FFT3(C->Nx, C->Ny, C->Nz, P -> P1_RowMaj, P -> P1_3Dk_RowMaj);
	FFT3(C->Nx, C->Ny, C->Nz, P -> P2_RowMaj, P -> P2_3Dk_RowMaj);
	FFT3(C->Nx, C->Ny, C->Nz, P -> P3_RowMaj, P -> P3_3Dk_RowMaj);
 */
	// --- Particular Solution --- //
	ParticularSolution(C, P, E, 0, C -> total);
	
	// k1 = k2 = k3 = 0 point, k-space origin
	E -> potential_A_3Dk_RowMaj[0] = 0; 
	E -> E_A_1_3Dk_RowMaj[0] = 0;
	E -> E_A_2_3Dk_RowMaj[0] = 0;
	E -> E_A_3_3Dk_RowMaj[0] = 0;

	// IFFT3 potential_A_3Dk_RowMaj
	IFFT3(C -> Nx, C -> Ny, C -> Nz, E -> potential_A_3Dk_RowMaj, E -> potential_A_RowMaj);
	
	// IFFT3 the 3Dk E Fields
	IFFT3(C -> Nx, C -> Ny, C -> Nz, E -> E_A_1_3Dk_RowMaj, E -> E_A_1_RowMaj);
	IFFT3(C -> Nx, C -> Ny, C -> Nz, E -> E_A_2_3Dk_RowMaj, E -> E_A_2_RowMaj);
	IFFT3(C -> Nx, C -> Ny, C -> Nz, E -> E_A_3_3Dk_RowMaj, E -> E_A_3_RowMaj);

	
	
	// --- Homo Solution --- //
	
	// BC @ interface_index + film_index
	for(int i = 0; i < C -> Nx; i++){
	for(int j = 0; j < C -> Ny; j++){
		E -> potential_bc_film_RowMaj     [j + C -> Ny * i] = 
			-(E -> potential_A_RowMaj[C -> film_index      + C -> Nz * (j + C -> Ny * i)]);
		E -> potential_bc_interface_RowMaj[j + C -> Ny * i] = 
			-(E -> potential_A_RowMaj[C -> interface_index + C -> Nz * (j + C -> Ny * i)]);
	}}
	
	// FFT2 the film and interface BC
	FFT2(C -> Nx, C -> Ny, E -> potential_bc_film_RowMaj, E -> potential_bc_film_2Dk_RowMaj);
	FFT2(C -> Nx, C -> Ny, E -> potential_bc_interface_RowMaj, E -> potential_bc_interface_2Dk_RowMaj);
	
	HomoSolution(C, E, S, 0, C -> Nx);
	
	// Kspace origin
	double complex bc[2] = {0, 0}; 
	double complex coeffs[2] = {0, 0};
	double complex one = 1.0; 
	double complex zero = 0.0;
	char N = 'N';
	int n2 = 2;
	int n1 = 1;
	
	// RHS [interface; film]
	bc[0] = E -> potential_bc_interface_2Dk_RowMaj[0];
	bc[1] = E -> potential_bc_film_2Dk_RowMaj[0];
	
	// bc_mat_inv * RHS 
	zgemm_(&N, &N, &n2, &n1, &n2, 
			&one, S -> electric_bc_mat_inv_korigin_ColMaj, &n2, 
			bc, &n2, 
			&zero, coeffs, &n2);
			
 	// copy the k1 = k2 = 0 point
	for(int z = 0; z < C -> Nz; z++)
	{
		E -> potential_B_2Dk_RowMaj[z] = coeffs[0] * C -> z_axis[z] + coeffs[1];
		E -> potential_B_d3_2Dk_RowMaj[z] = coeffs[0];
	}
	
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, 
		E -> potential_B_2Dk_RowMaj, E -> potential_B_RowMaj, 0, C -> Nz);
	
	// Find Potential = Potential_A + Potential_B
	for(int i = 0; i < C -> total; i++)
	{
		 T -> Potential_RowMaj[i] = E -> potential_B_RowMaj[i] + E -> potential_A_RowMaj[i]; 
	}	

	// Find E field
	for(int i = 0; i < C -> total; i++)
	{
		E -> E_B_1_2Dk_RowMaj[i] = -1 * I * C -> kx_grid_3D_RowMaj[i] * E -> potential_B_2Dk_RowMaj[i];
		E -> E_B_2_2Dk_RowMaj[i] = -1 * I * C -> ky_grid_3D_RowMaj[i] * E -> potential_B_2Dk_RowMaj[i];
		E -> E_B_3_2Dk_RowMaj[i] = -(E -> potential_B_d3_2Dk_RowMaj[i]);	
	}
	
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, E -> E_B_1_2Dk_RowMaj, E -> E_B_1_RowMaj, 0, C -> Nz);
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, E -> E_B_2_2Dk_RowMaj, E -> E_B_2_RowMaj, 0, C -> Nz);
	IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, E -> E_B_3_2Dk_RowMaj, E -> E_B_3_RowMaj, 0, C -> Nz);	

	for(int i = 0; i < C -> total; i++)
	{
		T -> E_1_depol_RowMaj[i] = E -> E_A_1_RowMaj[i] + E -> E_B_1_RowMaj[i];
		T -> E_2_depol_RowMaj[i] = E -> E_A_2_RowMaj[i] + E -> E_B_2_RowMaj[i];
		T -> E_3_depol_RowMaj[i] = E -> E_A_3_RowMaj[i] + E -> E_B_3_RowMaj[i];

		T -> f1_elec_RowMaj[i] = -0.5 * T -> E_1_depol_RowMaj[i];
		T -> f2_elec_RowMaj[i] = -0.5 * T -> E_2_depol_RowMaj[i];
		T -> f3_elec_RowMaj[i] = -0.5 * T -> E_3_depol_RowMaj[i];
	}
}

// Multithreaded code useless???

/*
	// Mulithreaded Version
	if(C -> MULTITHREAD)
	{
		ParticularSolution_MultiThreaded(C, P, E);
	}
	{
		ParticularSolution(C, P, E, 0, C -> total);
	}
	
	
	if(C -> MULTITHREAD)
	{
		HomoSolution_MultiThreaded(C, E, S);
	}
	else
	{
		HomoSolution(C, E, S, 0, C -> Nx);
	}
	
	
 	// IFFT2 the slices
	if(C -> MULTITHREAD)
	{
		IFFT2_slices_MultiThread(C -> Nx, C -> Ny, C -> Nz, 
			E -> potential_B_2Dk_RowMaj, E -> potential_B_RowMaj, C -> NUM_THREADS);		
	}
	else
	{
		IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, 
			E -> potential_B_2Dk_RowMaj, E -> potential_B_RowMaj, 0, C -> Nz);
	} 
	

 	// IFFT2
	if(C -> MULTITHREAD)
	{
		IFFT2_slices_MultiThread(C -> Nx, C -> Ny, C -> Nz, E -> E_B_1_2Dk_RowMaj, E -> E_B_1_RowMaj, C -> NUM_THREADS);
		IFFT2_slices_MultiThread(C -> Nx, C -> Ny, C -> Nz, E -> E_B_2_2Dk_RowMaj, E -> E_B_2_RowMaj, C -> NUM_THREADS);
		IFFT2_slices_MultiThread(C -> Nx, C -> Ny, C -> Nz, E -> E_B_3_2Dk_RowMaj, E -> E_B_3_RowMaj, C -> NUM_THREADS);		
	}
	else
	{
		IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, E -> E_B_1_2Dk_RowMaj, E -> E_B_1_RowMaj, 0, C -> Nz);
		IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, E -> E_B_2_2Dk_RowMaj, E -> E_B_2_RowMaj, 0, C -> Nz);
		IFFT2_slices(C -> Nx, C -> Ny, C -> Nz, E -> E_B_3_2Dk_RowMaj, E -> E_B_3_RowMaj, 0, C -> Nz);	
	} 

// Mulithreaded Version --> FASTER :)
typedef struct
{
	
	struct ConstantsStruct * C;
	struct PVector * P; 
	struct ElectricWorkspace * E;
	
	int thread_num;
	
} ParticularSolution_MultiThreadedStruct;

// Helper function that each thread runs
// Goes from f * i -> f * i - 1 in vector calculation
void * ParticularSolution_MultiThreaded_Helper(void * arg)
{
	ParticularSolution_MultiThreadedStruct * actual_arg = 
		(ParticularSolution_MultiThreadedStruct *) arg;
	struct ConstantsStruct * C = actual_arg -> C;
	struct PVector * P = actual_arg -> P;
	struct ElectricWorkspace * E = actual_arg -> E;

	int thread_num = actual_arg -> thread_num;

	// Where to calculate?
	int NUM_THREADS = C -> NUM_THREADS;
	int total = (C -> Nx) * (C -> Ny) * (C -> Nz);
	int f = (int) total / NUM_THREADS;
	int start = thread_num * f;
	int end = (thread_num + 1) * f;
	if(end >= total)
		end = total;

	ParticularSolution(C, P, E, start, end);	
}

// Main driver
void ParticularSolution_MultiThreaded(struct ConstantsStruct * C, 
	struct PVector * P, struct ElectricWorkspace * E)
{
	int NUM_THREADS = C -> NUM_THREADS;
	
	pthread_t thread[NUM_THREADS];
	int tid[NUM_THREADS];
	ParticularSolution_MultiThreadedStruct args[NUM_THREADS];
	
	for(int i = 0; i < NUM_THREADS; i++)
	{
		tid[i] = i;
		
		// argument construction
		args[i].C = C;
		args[i].P = P;
		args[i].E = E;
		args[i].thread_num = i;

		// create thread
		pthread_create(&thread[i], NULL, ParticularSolution_MultiThreaded_Helper, &args[i]);
	}
	
	// wait for threads to finish
	for(int i = 0; i < NUM_THREADS; i++)
		pthread_join(thread[i], NULL);
}


// Multithreaded
typedef struct
{
	
	struct ConstantsStruct * C;
	struct ElectricWorkspace * E;
	struct ElectricSolutionStruct_ColMaj * S;
	
	int thread_num;
	
} HomoSolution_MultiThreadedStruct;


void * HomoSolution_MultiThreaded_Helper(void * arg)
{
	HomoSolution_MultiThreadedStruct * actual_arg = 
		(HomoSolution_MultiThreadedStruct *) arg;
	struct ConstantsStruct * C = actual_arg -> C;
	struct ElectricWorkspace * E = actual_arg -> E;
	struct ElectricSolutionStruct_ColMaj * S = actual_arg -> S;

	int thread_num = actual_arg -> thread_num;

	// Where to calculate?
	int NUM_THREADS = C -> NUM_THREADS;
	int total = (C -> Nx);
	int f = (int) total / NUM_THREADS;
	int start = thread_num * f;
	int end = (thread_num + 1) * f;
	if(end >= total)
		end = total;

	HomoSolution(C, E, S, start, end);
}

void HomoSolution_MultiThreaded(struct ConstantsStruct * C, struct ElectricWorkspace * E, 
	struct ElectricSolutionStruct_ColMaj * S)
{
	int NUM_THREADS = C -> NUM_THREADS;
	
	pthread_t thread[NUM_THREADS];
	int tid[NUM_THREADS];
	HomoSolution_MultiThreadedStruct args[NUM_THREADS];
	
	for(int i = 0; i < NUM_THREADS; i++)
	{
		tid[i] = i;
		
		// argument construction
		args[i].C = C;
		args[i].E = E;
		args[i].S = S;
		args[i].thread_num = i;
	
		// create thread
		pthread_create(&thread[i], NULL, HomoSolution_MultiThreaded_Helper, &args[i]);
	}
	
	// wait for threads to finish
	for(int i = 0; i < NUM_THREADS; i++)
		pthread_join(thread[i], NULL);
}
*/