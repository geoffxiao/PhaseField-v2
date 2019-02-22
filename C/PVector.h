#ifndef PVECTOR_H
#define PVECTOR_H

#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "HelperFunctions.h"

struct PVector
{
	double complex * P1_RowMaj;
	double complex * P2_RowMaj;
	double complex * P3_RowMaj;

	double complex * P1_3Dk_RowMaj;
	double complex * P2_3Dk_RowMaj;
	double complex * P3_3Dk_RowMaj;
	
} PVector;

void FreePVector(struct PVector * P);
void InitPVector(int Nx, int Ny, int Nz, int type, struct PVector * P);

#endif