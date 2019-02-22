#ifndef ELECTRIC_SETUP_H
#define ELECTRIC_SETUP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Setup.h"
#include "HelperFunctions.h"

struct ElectricSolutionStruct{
	
	// double **** electric_bc_mat;
	double complex *** electric_bc_mat_inv_ColMaj;
	double complex * electric_bc_mat_inv_korigin_ColMaj;
	double ** electric_p_mat;
	
} ElectricSolutionStruct;

void ElectricSetup( struct ConstantsStruct * Constants,
					struct ElectricSolutionStruct * ElectricSolution);
					
void FreeElectricSolution( int Nx, int Ny, struct ElectricSolutionStruct * ElectricSolution);

#endif