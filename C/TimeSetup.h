#ifndef TIME_SETUP_H
#define TIME_SETUP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Setup.h"
#include "HelperFunctions.h"
#include "PVector.h"

struct TimeSolutionStruct{
	
	// First Order
	double complex ** time_q_1_1storder;
	double complex ** time_q_2_1storder;
	double complex ** time_q_3_1storder;
	
	double complex *** time_bc_mat_inv_1_1storder_ColMaj;
	double complex *** time_bc_mat_inv_2_1storder_ColMaj;
	double complex *** time_bc_mat_inv_3_1storder_ColMaj;
	
	// double ** time_scaling_1_1storder;
	// double ** time_scaling_2_1storder;
	// double ** time_scaling_3_1storder;
	
	// Second Order
	double complex ** time_q_1_2ndorder;
	double complex ** time_q_2_2ndorder;
	double complex ** time_q_3_2ndorder;
	
	double complex *** time_bc_mat_inv_1_2ndorder_ColMaj;
	double complex *** time_bc_mat_inv_2_2ndorder_ColMaj;
	double complex *** time_bc_mat_inv_3_2ndorder_ColMaj;
	
	// double ** time_scaling_1_2ndorder;
	// double ** time_scaling_2_2ndorder;
	// double ** time_scaling_3_2ndorder;
		
} TimeSolutionStruct;

void TimeSetup( struct ConstantsStruct * Constants, 
				struct TimeSolutionStruct * TimeSolution );

void FreeTimeSolution( int Nx, int Ny, struct TimeSolutionStruct * TimeSolution );

#endif