#include "CalcLandauEnergy.h"

void CalcLandauEnergy(struct ConstantsStruct * C, struct TotalStateStruct * T,
	struct PVector * P, int start, int end)
{
	double a1 = C->a1;
	double a11 = C->a11;
	double a12 = C->a12;
	double a111 = C->a111;
	double a112 = C->a112;
	double a123 = C->a123;
	double a1111 = C->a1111;
	double a1112 = C->a1112;
	double a1122 = C->a1122;
	double a1123 = C->a1123;
	
	double complex * P1 = P->P1_RowMaj;
	double complex * P2 = P->P2_RowMaj;
	double complex * P3 = P->P3_RowMaj;
	
	double complex P1_val, P2_val, P3_val;
	for(int i = start; i < end; i++)
	{
		P1_val = P1[i];
		P2_val = P2[i];
		P3_val = P3[i];
		
		T->f1_landau_RowMaj[i] = ( P1_val * (
			2.0 * a1 + 4.0 * a11 * P1_val * P1_val + 2.0 * a12 * ( P2_val * P2_val + P3_val * P3_val ) +
			6.0 * a111 * P1_val * P1_val * P1_val * P1_val + 2.0 * a112 * ( P2_val * P2_val * P2_val * P2_val + P3_val * P3_val * P3_val * P3_val + 2.0 * P1_val * P1_val * ( P2_val * P2_val + P3_val * P3_val ) ) +
			2.0 * a123 * P2_val * P2_val * P3_val * P3_val +
			8.0 * a1111 * P1_val * P1_val * P1_val * P1_val * P1_val * P1_val + 2.0 * a1112* ( P2_val * P2_val * P2_val * P2_val * P2_val * P2_val + P3_val * P3_val * P3_val * P3_val * P3_val * P3_val + 3 * P1_val * P1_val * P1_val * P1_val * ( P2_val * P2_val + P3_val * P3_val ) ) +
			4.0 * a1122 * P1_val * P1_val * ( P2_val * P2_val * P2_val * P2_val + P3_val * P3_val * P3_val * P3_val ) +
			2.0 * a1123 * ( 2.0 * P1_val * P1_val * P2_val * P2_val * P3_val * P3_val + P2_val * P2_val * P3_val * P3_val * ( P2_val * P2_val + P3_val * P3_val ) ) ) );

		T->f2_landau_RowMaj[i] = ( P2_val * (
			2.0 * a1 + 4.0 * a11 * P2_val * P2_val + 2.0 * a12 * ( P3_val * P3_val + P1_val * P1_val ) +
			6.0 * a111 * P2_val * P2_val * P2_val * P2_val + 2.0 * a112 * ( P3_val * P3_val * P3_val * P3_val + P1_val * P1_val * P1_val * P1_val + 2.0 * P2_val * P2_val * ( P3_val * P3_val + P1_val * P1_val ) ) +
			2.0 * a123 * P3_val * P3_val * P1_val * P1_val +
			8.0 * a1111 * P2_val * P2_val * P2_val * P2_val * P2_val * P2_val + 2.0 * a1112* ( P3_val * P3_val * P3_val * P3_val * P3_val * P3_val + P1_val * P1_val * P1_val * P1_val * P1_val * P1_val + 3 * P2_val * P2_val * P2_val * P2_val * ( P3_val * P3_val + P1_val * P1_val ) ) +
			4.0 * a1122 * P2_val * P2_val * ( P3_val * P3_val * P3_val * P3_val + P1_val * P1_val * P1_val * P1_val ) +
			2.0 * a1123 * ( 2.0 * P2_val * P2_val * P3_val * P3_val * P1_val * P1_val + P3_val * P3_val * P1_val * P1_val * ( P3_val * P3_val + P1_val * P1_val ) ) ) );

		T->f3_landau_RowMaj[i] = ( P3_val * (
			2.0 * a1 + 4.0 * a11 * P3_val * P3_val + 2.0 * a12 * ( P1_val * P1_val + P2_val * P2_val ) +
			6.0 * a111 * P3_val * P3_val * P3_val * P3_val + 2.0 * a112 * ( P1_val * P1_val * P1_val * P1_val + P2_val * P2_val * P2_val * P2_val + 2.0 * P3_val * P3_val * ( P1_val * P1_val + P2_val * P2_val ) ) +
			2.0 * a123 * P1_val * P1_val * P2_val * P2_val +
			8.0 * a1111 * P3_val * P3_val * P3_val * P3_val * P3_val * P3_val + 2.0 * a1112* ( P1_val * P1_val * P1_val * P1_val * P1_val * P1_val + P2_val * P2_val * P2_val * P2_val * P2_val * P2_val + 3 * P3_val * P3_val * P3_val * P3_val * ( P1_val * P1_val + P2_val * P2_val ) ) +
			4.0 * a1122 * P3_val * P3_val * ( P1_val * P1_val * P1_val * P1_val + P2_val * P2_val * P2_val * P2_val ) +
			2.0 * a1123 * ( 2.0 * P3_val * P3_val * P1_val * P1_val * P2_val * P2_val + P1_val * P1_val * P2_val * P2_val * ( P1_val * P1_val + P2_val * P2_val ) ) ) );
	}
}



// Mulithreading...
typedef struct
{
	
	struct ConstantsStruct * C;
	struct PVector * P; 
	struct TotalStateStruct * T;
	
	int thread_num;
	
} CalcLandauEnergy_MultiThreadedStruct;

// Helper function that each thread runs
// Goes from f * i -> f * i - 1 in vector calculation
void * CalcLandauEnergy_MultiThreaded_Helper(void * arg)
{
	CalcLandauEnergy_MultiThreadedStruct * actual_arg = 
		(CalcLandauEnergy_MultiThreadedStruct *) arg;
	struct ConstantsStruct * C = actual_arg -> C;
	struct PVector * P = actual_arg -> P;
	struct TotalStateStruct * T = actual_arg -> T;

	int thread_num = actual_arg -> thread_num;

	// Where to calculate?
	int NUM_THREADS = C -> NUM_THREADS;
	int total = (C -> Nx) * (C -> Ny) * (C -> Nz);
	int f = (int) total / NUM_THREADS;
	int start = thread_num * f;
	int end = (thread_num + 1) * f;
	if(end >= total)
		end = total;

	CalcLandauEnergy(C, T, P, start, end);	
}

// Main driver
void CalcLandauEnergy_MultiThreaded(struct ConstantsStruct * C, struct TotalStateStruct * T,
	struct PVector * P)
{
	int NUM_THREADS = C -> NUM_THREADS;
	
	pthread_t thread[NUM_THREADS];
	int tid[NUM_THREADS];
	CalcLandauEnergy_MultiThreadedStruct args[NUM_THREADS];
	
	for(int i = 0; i < NUM_THREADS; i++)
	{
		tid[i] = i;
		
		// argument construction
		args[i].C = C;
		args[i].P = P;
		args[i].T = T;
		args[i].thread_num = i;

		// create thread
		pthread_create(&thread[i], NULL, CalcLandauEnergy_MultiThreaded_Helper, &args[i]);
	}
	
	// wait for threads to finish
	for(int i = 0; i < NUM_THREADS; i++)
		pthread_join(thread[i], NULL);
}