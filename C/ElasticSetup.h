#ifndef ELASTIC_SETUP_H
#define ELASTIC_SETUP_H

#include <stdio.h>
#include <stdlib.h>

#include "Setup.h"
#include "HelperFunctions.h"

#ifndef M_PI
	#define M_PI (3.14159265358979323846264338327950288)
#endif

struct ElasticSolutionStruct{
	
	double **** elastic_eigenproblem_mat;
	
	// double complex *** elastic_scaling_mat;
	double complex *** elastic_bc_mat_inv_ColMaj;
	double complex * elastic_bc_mat_inv_korigin_ColMaj;
	double complex **** elastic_eigenvec_mat; 
	double complex *** elastic_eigenval_mat;
	double complex *** elastic_bc_mat_L_ColMaj;
	double complex *** elastic_bc_mat_L_inv_ColMaj;
	double complex *** elastic_bc_mat_U_ColMaj;
	double complex *** elastic_bc_mat_P_ColMaj;
	double complex **** elastic_bc_mat;
	
	double complex * Green_11_RowMaj;
	double complex * Green_12_RowMaj;
	double complex * Green_13_RowMaj;
	
	double complex * Green_21_RowMaj;
	double complex * Green_22_RowMaj;
	double complex * Green_23_RowMaj;
	
	double complex * Green_31_RowMaj;
	double complex * Green_32_RowMaj;
	double complex * Green_33_RowMaj;
	
} ElasticSolutionStruct;

void ElasticSetup( struct ConstantsStruct * Constants,
					struct ElasticSolutionStruct * ElasticSolution );
					
void FreeElasticSolution( int Nx, int Ny, 
	struct ElasticSolutionStruct * ElasticSolution );
	
#endif