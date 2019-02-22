#ifndef ENERGIES_H
#define ENERGIES_H

#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct EnergiesStruct
{
	double complex * f1_RowMaj;
	double complex * f2_RowMaj;
	double complex * f3_RowMaj;

	double complex * f1_3Dk_RowMaj;
	double complex * f2_3Dk_RowMaj;
	double complex * f3_3Dk_RowMaj;
	
} EnergiesStruct;

void FreeEnergies(struct EnergiesStruct * F);
void InitEnergies(int Nx, int Ny, int Nz, struct EnergiesStruct * F);

#endif