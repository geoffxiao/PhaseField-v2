#ifndef ELECTRICWORKSPACE_H
#define ELECTRICWORKSPACE_H

#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct ElectricWorkspaceStruct
{

	double complex * potential_A_3Dk_RowMaj;
	double complex * potential_A_RowMaj;
	
	double complex * potential_B_2Dk_RowMaj;
	double complex * potential_B_RowMaj;	
	double complex * potential_B_d3_2Dk_RowMaj;
	
	double complex * E_A_1_3Dk_RowMaj;
	double complex * E_A_2_3Dk_RowMaj;
	double complex * E_A_3_3Dk_RowMaj;
	
	double complex * E_A_1_RowMaj;
	double complex * E_A_2_RowMaj;
	double complex * E_A_3_RowMaj;
	
	double complex * E_B_1_2Dk_RowMaj;
	double complex * E_B_2_2Dk_RowMaj;
	double complex * E_B_3_2Dk_RowMaj;
	
	double complex * E_B_1_RowMaj;
	double complex * E_B_2_RowMaj;
	double complex * E_B_3_RowMaj;
	
	double complex * potential_bc_interface_RowMaj;
	double complex * potential_bc_film_RowMaj;
	
	double complex * potential_bc_interface_2Dk_RowMaj;
	double complex * potential_bc_film_2Dk_RowMaj;	

} ElectricWorkspaceStruct;


void FreeElectricWorkspace(struct ElectricWorkspaceStruct * E);
void InitElectricWorkspace(int Nx, int Ny, int Nz, struct ElectricWorkspaceStruct * E);

#endif