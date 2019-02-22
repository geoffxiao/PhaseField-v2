#ifndef EIGENSTRAIN_H
#define EIGENSTRAIN_H

#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "HelperFunctions.h"
#include "Setup.h"
#include "PVector.h"

struct EigenstrainStruct
{
	
	double complex * Eigenstrain_11_RowMaj;
	double complex * Eigenstrain_22_RowMaj;
	double complex * Eigenstrain_33_RowMaj;
	double complex * Eigenstrain_12_RowMaj;	
	double complex * Eigenstrain_13_RowMaj;
	double complex * Eigenstrain_23_RowMaj;
	
	double complex * Eigenstrain_11_3Dk_RowMaj;
	double complex * Eigenstrain_22_3Dk_RowMaj;
	double complex * Eigenstrain_33_3Dk_RowMaj;
	double complex * Eigenstrain_12_3Dk_RowMaj;	
	double complex * Eigenstrain_13_3Dk_RowMaj;
	double complex * Eigenstrain_23_3Dk_RowMaj;	
	
} EigenstrainStruct;

void InitEigenstrain(int Nx, int Ny, int Nz, struct EigenstrainStruct * Eigenstrain);
void CalcEigenstrain(struct EigenstrainStruct * Eigenstrain, 
	struct PVector * P, struct ConstantsStruct * Constants);
void FreeEigenstrain(struct EigenstrainStruct * Eigenstrain);

#endif