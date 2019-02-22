#ifndef TOTAL_STATE_H
#define TOTAL_STATE_H

#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct TotalStateStruct
{
	
	double complex * Potential_RowMaj;
	
	double complex * E_1_depol_RowMaj;
	double complex * E_2_depol_RowMaj;
	double complex * E_3_depol_RowMaj;
	
	double complex * f1_elec_RowMaj;
	double complex * f2_elec_RowMaj;
	double complex * f3_elec_RowMaj;
	
	double complex * f1_elastic_RowMaj;
	double complex * f2_elastic_RowMaj;
	double complex * f3_elastic_RowMaj;
	
	double complex * u_1_RowMaj;
	double complex * u_2_RowMaj;
	double complex * u_3_RowMaj;
	
	double complex * e_11_RowMaj;
	double complex * e_12_RowMaj;
	double complex * e_13_RowMaj;
	double complex * e_22_RowMaj;
	double complex * e_23_RowMaj;
	double complex * e_33_RowMaj;

	double complex * TotalStrain_11_RowMaj;
	double complex * TotalStrain_12_RowMaj;
	double complex * TotalStrain_13_RowMaj;
	double complex * TotalStrain_22_RowMaj;
	double complex * TotalStrain_23_RowMaj;
	double complex * TotalStrain_33_RowMaj;
	
	double complex * f1_landau_RowMaj;
	double complex * f2_landau_RowMaj;
	double complex * f3_landau_RowMaj;
	
	
} TotalStateStruct;

void FreeTotalState(struct TotalStateStruct * TotalState);
void InitTotalState(int Nx, int Ny, int Nz, struct TotalStateStruct * TotalState);

#endif