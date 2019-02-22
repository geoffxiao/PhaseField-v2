#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "HelperFunctions.h"
#include "Setup.h"
#include "ElasticSetup.h"
#include "ElectricSetup.h"
#include "TimeSetup.h"
#include "CalcElecEnergy.h"
#include "CalcElasticEnergy.h"
#include "PVector.h"
#include "ElectricWorkspace.h"
#include "TotalState.h"
#include "Eigenstrain.h"
#include "ElasticWorkspace.h"
#include "Energies.h"
#include "TimeStep.h"
#include "CalcLandauEnergy.h"


// LAPACK uses column major matrices
// zgemm, dgemm, etc...
// FFTW uses row major
// row major = transpose(column major)
// No memory allocation done in Main
// All dynamic memory alloc done in the helpers

// FFTW uses row major matrices

// Main function
int main(int argc, char ** argv)
{
	// Set up randomization
	srand(1);


	
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	// --- Set up Constants --- //
	struct ConstantsStruct Constants;
	InitSetup(&Constants);

	int Nx = Constants.Nx;
	int Ny = Constants.Ny;
	int Nz = Constants.Nz;
	
	int total = Constants.total;
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //



	
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	// --- Elastic Setup --- //
	struct ElasticSolutionStruct ElasticSolution;
   	ElasticSetup(&Constants, &ElasticSolution);
	// elastic_eigenvec_mat[k1][k2][values in the eigenvec][eigenvector #]
	// 		Same as MATLAB code ordering!
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //

	
 	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //	
	// --- Electric Setup --- //
	struct ElectricSolutionStruct ElectricSolution;
	ElectricSetup(&Constants, &ElectricSolution);
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	
	
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	// --- Time Setup --- //
	struct TimeSolutionStruct TimeSolution;
	TimeSetup(&Constants, &TimeSolution);
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //	
	
	

	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	// Other misc. init
	
	// Workspaces to avoid repeated malloc and free... is it necessary?
	struct ElectricWorkspaceStruct ElectricWorkspace;
	InitElectricWorkspace(Nx, Ny, Nz, &ElectricWorkspace);
	
	struct ElasticWorkspace U;
	InitElasticWorkspace(Nx, Ny, Nz, &U);

	// Store the total state
	// f elec, f elastic, f landau
	// u_i, totalstrain_ij
	// E_depol
	// potential V
	struct TotalStateStruct TotalState;
	InitTotalState(Nx, Ny, Nz, &TotalState);
	
	// eigenstrains
	struct EigenstrainStruct Eigenstrain;
	InitEigenstrain(Nx, Ny, Nz, &Eigenstrain);
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //


	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	// Init
	// P(t = 0)
	// 0 = randomize
	// 1 = load
	struct PVector P_t0; 
	InitPVector(Nx, Ny, Nz, 0, &P_t0);
	
	// F(t = 0), total energy at t = 0
	struct EnergiesStruct F_t0;
	InitEnergies(Nx, Ny, Nz, &F_t0);
	
	// P(t = 1), the next P time step
	// 2 = zero P vector completely
	struct PVector P_t1;
	InitPVector(Nx, Ny, Nz, 2, &P_t1);
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //

	
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	// Main Driver
	int c = 0; int save_c = 0;
	char c_char[256];
	double Err1, Err2, Err3, Err;
	Err = 1e10;

	while(c < Constants.MaxIter && Err > Constants.epsilon)
	{
		// Find Energies for P_t0
		// electrical, elastic, and landau
		// General notes:
		//		P_t0 includes FFT3'ed and real space P already
		//		Store impt stuff into TotalState struct
		CalcElecEnergy(&P_t0, &Constants, &ElectricSolution, &ElectricWorkspace, &TotalState);
		CalcElasticEnergy(&U, &P_t0, &Eigenstrain, &Constants, &ElasticSolution, &TotalState);
		CalcLandauEnergy_MultiThreaded(&Constants, &TotalState, &P_t0);
		
		// Find Total Energy
		for(int i = 0; i < total; i++)
		{
			F_t0.f1_RowMaj[i] = (TotalState.f1_elec_RowMaj[i] + TotalState.f1_elastic_RowMaj[i] + TotalState.f1_landau_RowMaj[i]) * Constants.in_film_RowMaj[i];
			F_t0.f2_RowMaj[i] = (TotalState.f2_elec_RowMaj[i] + TotalState.f2_elastic_RowMaj[i] + TotalState.f2_landau_RowMaj[i]) * Constants.in_film_RowMaj[i];
			F_t0.f3_RowMaj[i] = (TotalState.f3_elec_RowMaj[i] + TotalState.f3_elastic_RowMaj[i] + TotalState.f3_landau_RowMaj[i]) * Constants.in_film_RowMaj[i];
		}
		// FFT3 the energies
		FFT3(Nx, Ny, Nz, F_t0.f1_RowMaj, F_t0.f1_3Dk_RowMaj);
		FFT3(Nx, Ny, Nz, F_t0.f2_RowMaj, F_t0.f2_3Dk_RowMaj);
		FFT3(Nx, Ny, Nz, F_t0.f3_RowMaj, F_t0.f3_3Dk_RowMaj);

		// Find P_t1 from P_t0
		//		Also FFT3 the calculated P
		TimeStep_1stOrd(&P_t0, &F_t0, &Constants, &TimeSolution, &P_t1);
		
		// Calculate Error
		Err1 = CalcError(P_t0.P1_RowMaj, P_t1.P1_RowMaj, total);
		Err2 = CalcError(P_t0.P2_RowMaj, P_t1.P2_RowMaj, total);
		Err3 = CalcError(P_t0.P3_RowMaj, P_t1.P3_RowMaj, total);
		Err = Err1 + Err2 + Err3;
		
		// Let P_t0 = P_t1
		free(P_t0.P1_RowMaj); free(P_t0.P1_3Dk_RowMaj);
		free(P_t0.P2_RowMaj); free(P_t0.P2_3Dk_RowMaj);
		free(P_t0.P3_RowMaj); free(P_t0.P3_3Dk_RowMaj);

		P_t0.P1_RowMaj = P_t1.P1_RowMaj;
		P_t0.P2_RowMaj = P_t1.P2_RowMaj;
		P_t0.P3_RowMaj = P_t1.P3_RowMaj;

		P_t0.P1_3Dk_RowMaj = P_t1.P1_3Dk_RowMaj;
		P_t0.P2_3Dk_RowMaj = P_t1.P2_3Dk_RowMaj;
		P_t0.P3_3Dk_RowMaj = P_t1.P3_3Dk_RowMaj;
						

		// Rezero P_t1
		InitPVector(Nx, Ny, Nz, 2, &P_t1);
		
		if(c == Constants.saves[save_c])
		{
			// Save into directory
			sprintf(c_char, "./%d", c);
			mkdir(c_char, 0700);
			
			sprintf(c_char, "./%d/P1_%d", c, c);
			SaveFile_Complex(c_char, total, 0, P_t0.P1_RowMaj);
			
			sprintf(c_char, "./%d/P2_%d", c, c);
			SaveFile_Complex(c_char, total, 0, P_t0.P2_RowMaj);
					
			sprintf(c_char, "./%d/P3_%d", c, c);
			SaveFile_Complex(c_char, total, 0, P_t0.P3_RowMaj);
			
			save_c++;
		}
		
		// Progress
		printf("Step: %d\nError: %.4e\n\n", c, Err);
		c++;
	}
	
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //
	// Save the final outputs
	mkdir("./final", 0700);
	SaveFile_Complex("./final/P1_final.txt", total, 0, P_t0.P1_RowMaj);
	SaveFile_Complex("./final/P2_final.txt", total, 0, P_t0.P2_RowMaj);
	SaveFile_Complex("./final/P3_final.txt", total, 0, P_t0.P3_RowMaj);
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //	


	
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //	
 	// Memory clean up
	FreeElasticSolution(Nx, Ny, &ElasticSolution);
	FreeElectricSolution(Nx, Ny, &ElectricSolution); 	
	FreeTimeSolution(Nx, Ny, &TimeSolution);
	
	// Free P Vectors
	FreePVector(&P_t0); FreePVector(&P_t1);

	// Free Energy
	FreeEnergies(&F_t0);
	
	// Free Workspaces
	FreeElectricWorkspace(&ElectricWorkspace);
	FreeElasticWorkspace(&U);
	
	// Free structs
	FreeSetup(&Constants);
	FreeTotalState(&TotalState);
	FreeEigenstrain(&Eigenstrain);
	// --------------------------------------------------------------------- //
	// --------------------------------------------------------------------- //	

	
	return 0;
}