// Make sure only defined once
#ifndef SETUP_H
#define SETUP_H

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "HelperFunctions.h"

typedef int bool;
#define true 1
#define false 0

#ifndef M_PI
	#define M_PI (3.14159265358979323846264338327950288)
#endif

// Constants
struct ConstantsStruct{

	double Temperature;
	double dt_factor;
	double ThermConst;
	double epsilon;
	double BTO_pct;

	int Nx, Ny, Nz;
	int total;
	int interface_index, sub_index, film_index;
	int Nz_film;

	double E_1_applied, E_2_applied, E_3_applied;
	double Us_11, Us_22, Us_12;

	bool LoadFromFile;
	char * Pathname;
	char * Filename;
	
	double C11, C12, C44;
	double **** C;
	double a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, a1123;
	double b11, b12; // renorm a11, a12 b/c strain effects
	double q11, q12, q44;
	double Q11, Q12, Q44;
	double Lx, Ly, Lz;	
	double G11, G12, H14, H44, G110;
	double k_electric_11, k_electric_22, k_electric_33;
	
	double dt;
	int * saves;
	int MaxIter;

	double permittivity_0;
	double boltzmann;
	
	double * x_axis, * y_axis, * z_axis;
	double * kx_axis, * ky_axis, * kz_axis;
	
	double * kx_grid_3D_RowMaj, * ky_grid_3D_RowMaj, * kz_grid_3D_RowMaj;
	double * kx_grid_2D_RowMaj, * ky_grid_2D_RowMaj;
	double ** kx_grid_2D, ** ky_grid_2D;
	
	double complex * in_film_RowMaj;
	
	double dx, dy, dz;
	
	double h_film, h_sub, h_int;
	
	int NUM_THREADS;
	bool MULTITHREAD;
	
} ConstantsStruct;

void Setup_k_axis(double L, int N, double * k);

void InitSetup( struct ConstantsStruct * Constants );
void FreeSetup( struct ConstantsStruct * Constants );
double coth( double x );

#endif