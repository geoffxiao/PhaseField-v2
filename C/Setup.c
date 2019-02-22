#include "Setup.h"

// Initialize the Constants
void InitSetup( struct ConstantsStruct * Constants )
{
	// doubles have to have .0
	double Temperature = 25.0;
	double dt_factor = 0.05;
	double BTO_pct = 1.0;
	double ThermConst = 0.0;
	double epsilon = 1e-5;
	int MaxIter = 100;
	int * saves = malloc(sizeof(int) * MaxIter);
	
	for(int i = 0; i < MaxIter; i++)
	{
		saves[i] = 100;
	}
	int NUM_THREADS = 8;
	bool MULTITHREAD = 0;
	
	int Nx = 32;
	int Ny = 32;
	int Nz = 32;
	
	// sub_index
	// interface_index
	// film_index
	
	int sub_index = 2 - 1; // -1 b/c C starts with 0, how many indices is the substrate?
	int interface_index = sub_index + 1;
	int film_index = Nz - 4 - 1; // -1 b/c C starts with 0, how thick is the air?
	int Nz_film = film_index - sub_index;
	
	bool LoadFromFile = false;
	
	double E_1_applied = 0;
	double E_2_applied = 0;
	double E_3_applied = 0;		
	
	// 0 = SSO
	// 1 = GSO
	// 2 = NSO
	// -1 = Set Us_11, Us_22, Us_12 manually
	int substrate = 0; 
	
	double Us_11, Us_12, Us_22;
	
	// --- Setup --- //
	double SSO_lattice_11 = 3.9912;
	double GSO_lattice_11 = 3.9701;
	double NSO_lattice_11 = 4.0143;

	double SSO_lattice_22 = 3.9814;
	double GSO_lattice_22 = 3.9657;
	double NSO_lattice_22 = 4.0017;

	// http://www.mmm.psu.edu/JWang2010_JAP_BTO1stPrinciple.pdf
	double BTO_lattice = 4.006;
	// https://deepblue.lib.umich.edu/bitstream/handle/2027.42/62658/nature02773.pdf?sequence=1&isAllowed=y
	double STO_lattice = 3.905;
	double STO_pct = 1 - BTO_pct;
	double BST_lattice = BTO_pct * BTO_lattice + STO_pct * STO_lattice;

	if(substrate != -1)
	{
		double Substrate_lattice_11;
		double Substrate_lattice_22;
		switch(substrate)
		{
			case 0 : 
				Substrate_lattice_11 = SSO_lattice_11;
				Substrate_lattice_22 = SSO_lattice_22;
				break;
			
			case 1 :
				Substrate_lattice_11 = GSO_lattice_11;
				Substrate_lattice_22 = GSO_lattice_22;
				break;
			
			case 2 :
				Substrate_lattice_11 = NSO_lattice_11;
				Substrate_lattice_22 = NSO_lattice_22;
				break;
		}
		// http://www.iue.tuwien.ac.at/phd/dhar/node12.html	
		Us_11 = ((Substrate_lattice_11 - BST_lattice) / BST_lattice); // unitless misfit strain
		Us_22 = ((Substrate_lattice_22 - BST_lattice) / BST_lattice); // anisotropic misfit strain
		Us_12 = 0;
	}
	
	char * Pathname = "/";
	char Filename[256];
	sprintf(Filename, "BTO %g-%g", BTO_pct*100, STO_pct*100);


	// --- BTO Constants --- //
	// BTO =  http://www.ems.psu.edu/~chen/publications/YL2005APL.pdf
	double C11_BTO = 1.78 * 1e11;
	double C12_BTO = 0.964 * 1e11;
	double C44_BTO = 1.22 * 1e11;

	double Q11_BTO = 0.10; 
	double Q12_BTO = -0.034;
	double Q44_BTO = 0.029;

	double a1_BTO = 4.124 * ( Temperature - 115.0 ) * 1e5;
	double a11_BTO = -2.097 * 1e8;
	double a12_BTO = 7.974 * 1e8;
	double a111_BTO = 1.294 * 1e9;
	double a112_BTO = -1.950 * 1e9;
	double a123_BTO = -2.5 * 1e9;
	double a1111_BTO = 3.863 * 1e10;
	double a1112_BTO = 2.529 * 1e10;
	double a1122_BTO = 1.637 * 1e10;
	double a1123_BTO = 1.367 * 1e10;

	// --- STO Constants --- //
	double C11_STO = 3.36 * 1e11;
	double C12_STO = 1.07 * 1e11;
	double C44_STO = 1.27 * 1e11;

	// STO Yu Luan Li PRB 184112, Table II
	// 1 cm^4 / esu^2 = 8.988 * 1e10 m^4 / C^2
	// http://www.wolframalpha.com/input/?i=(1+cm)%5E4%2F(3.3356e-10+coulomb)%5E2
	double Q11_STO = 4.57 * 1e-2; 
	double Q12_STO = -1.348 * 1e-2;
	double Q44_STO = 0.957 * 1e-2;

	// STO Yu Luan Li PRB 184112, Table III a1 (1st, from ref 6, 16, 22), Table
	// VII a11, a12, a1
	// Shirokov, Yuzyu, PRB 144118, Table I STO in parantheses

	// 1 cm^2 * dyn / esu^2 = 8.9876 * 1e9 J * m^2 / C
	double a1_STO = 4.05 * 1e7 * ( coth( 54.0 / (Temperature+273.0) ) - coth(54.0/30.0) );

	// 1 cm^6 * dyn / esu^4 = 8.0776 * 1e20 J * m^5 / C^4
	double a11_STO = 17 * 1e8;
	double a12_STO = 13.7 * 1e8;

	double a111_STO = 0;
	double a112_STO = 0;
	double a123_STO = 0;
	double a1111_STO = 0;
	double a1112_STO = 0;
	double a1122_STO = 0;
	double a1123_STO = 0;
	
	// --- Linear Combination of Constants --- //
	double C11 = C11_BTO * BTO_pct + C11_STO * STO_pct;
	double C12 = C12_BTO * BTO_pct + C12_STO * STO_pct;
	double C44 = C44_BTO * BTO_pct + C44_STO * STO_pct;

	double q11_BTO = C11_BTO * Q11_BTO + 2.0 * C12_BTO * Q12_BTO;
	double q12_BTO = C11_BTO * Q12_BTO + C12_BTO * (Q11_BTO + Q12_BTO);
	double q44_BTO = 2.0 * C44_BTO * Q44_BTO;

	double q11_STO = C11_STO * Q11_STO + 2.0 * C12_STO * Q12_STO;
	double q12_STO = C11_STO * Q12_STO + C12_STO * (Q11_STO + Q12_STO);
	double q44_STO = 2.0 * C44_STO * Q44_STO;

	double q11 = q11_BTO * BTO_pct + q11_STO * STO_pct;
	double q12 = q12_BTO * BTO_pct + q12_STO * STO_pct;
	double q44 = q44_BTO * BTO_pct + q44_STO * STO_pct;

	double q11_h = q11 + 2.0 * q12;
	double q22_h = q11 - q12;
	double C11_h = C11 + 2.0 * C12;
	double C22_h = C11 - C12;

	double Q11 = (1.0/3.0) * ( (q11_h/C11_h) + (2.0 * q22_h/C22_h) );
	double Q12 = (1.0/3.0) * ( (q11_h/C11_h) - (q22_h/C22_h) );
	double Q44 = q44 / (2.0 * C44);
	
	double a1 = a1_BTO * BTO_pct + a1_STO * STO_pct;
	double a11 = a11_BTO * BTO_pct + a11_STO * STO_pct;
	double a12 = a12_BTO * BTO_pct + a12_STO * STO_pct;
	double a111 = a111_BTO * BTO_pct + a111_STO * STO_pct;
	double a112 = a112_BTO * BTO_pct + a112_STO * STO_pct;
	double a123 = a123_BTO * BTO_pct + a123_STO * STO_pct;
	double a1111 = a1111_BTO * BTO_pct + a1111_STO * STO_pct;
	double a1112 = a1112_BTO * BTO_pct + a1112_STO * STO_pct;
	double a1122 = a1122_BTO * BTO_pct + a1122_STO * STO_pct;
	double a1123 = a1123_BTO * BTO_pct + a1123_STO * STO_pct;

	double l_0 = 1e-9; double l_z = 1e-9;
	double Lx = l_0*(Nx-1); double Ly = l_0*(Ny-1); double Lz = l_z*(Nz-1);

	double a0 = abs( 4.124 * ( Temperature - 115.0 ) * 1e5 );
	double G110 = l_0 * l_0 * a0;
	
	double G11 = 1.0 * G110;
	double G12 = 0;
	double H14 = 0;
	double H44 = G11;
	
	double dt = dt_factor / a0;
	
	double k_electric_11 = 100.0;
	double k_electric_22 = 100.0;
	double k_electric_33 = 100.0;
	
	// Set up the real space axes
	double dx = Lx / (Nx - 1);
	double dy = Ly / (Ny - 1);
	double dz = Lz / (Nz - 1);
	
	double * x_axis = malloc(sizeof(double) * Nx);
	double * y_axis = malloc(sizeof(double) * Ny);
	double * z_axis = malloc(sizeof(double) * Nz);
	
	for(int i = 0; i < Nx; i++){ x_axis[i] = i * dx; }
	for(int i = 0; i < Ny; i++){ y_axis[i] = i * dy; }
	for(int i = 0; i < Nz; i++){ z_axis[i] = i * dz; }
	
	// Rescale z axis
	double rescale = z_axis[(int)(round(Nz / 2) - 1)];
	for(int i = 0; i < Nz; i++){ z_axis[i] = z_axis[i] - rescale; }
	
	// k-space	
	double * kx_axis = malloc(sizeof(double) * Nx);
	double * ky_axis = malloc(sizeof(double) * Ny);
	double * kz_axis = malloc(sizeof(double) * Nz);
	
 	Setup_k_axis(Lx, Nx, kx_axis);
	Setup_k_axis(Ly, Ny, ky_axis);
	Setup_k_axis(Lz, Nz, kz_axis);

/* 	double *** x_grid = AllocLargeArray_3D(Nx, Ny, Nz);
 	double *** y_grid = AllocLargeArray_3D(Nx, Ny, Nz);
	double *** z_grid = AllocLargeArray_3D(Nx, Ny, Nz);

	Meshgrid(Nx, Ny, Nz, x_axis, y_axis, z_axis, x_grid, y_grid, z_grid);
	 */

	double *** kx_grid_3D = AllocLargeArray_3D(Nx, Ny, Nz);
 	double *** ky_grid_3D = AllocLargeArray_3D(Nx, Ny, Nz);
	double *** kz_grid_3D = AllocLargeArray_3D(Nx, Ny, Nz);
	Meshgrid(Nx, Ny, Nz, kx_axis, ky_axis, kz_axis, 
			kx_grid_3D, ky_grid_3D, kz_grid_3D); 

	// Convert the 3D and 2D grids into a row major 1D array
	double * kx_grid_3D_RowMaj = malloc(sizeof(double) * Nx * Ny * Nz);
	double * ky_grid_3D_RowMaj = malloc(sizeof(double) * Nx * Ny * Nz);
	double * kz_grid_3D_RowMaj = malloc(sizeof(double) * Nx * Ny * Nz);
	
	Convert3D_2_RowMajor(Nx, Ny, Nz, kx_grid_3D, kx_grid_3D_RowMaj);
	Convert3D_2_RowMajor(Nx, Ny, Nz, ky_grid_3D, ky_grid_3D_RowMaj);
	Convert3D_2_RowMajor(Nx, Ny, Nz, kz_grid_3D, kz_grid_3D_RowMaj);
	
 	// Deep copy here
	double ** kx_grid_2D = AllocLargeArray_2D(Nx, Ny);
	double ** ky_grid_2D = AllocLargeArray_2D(Nx, Ny);
	for(int i = 0; i < Nx; i++)
	for(int j = 0; j < Ny; j++)
	{
		kx_grid_2D[i][j] = kx_grid_3D[i][j][0];
		ky_grid_2D[i][j] = ky_grid_3D[i][j][0];
	}

	Free_3D_mat(Nx, Ny, Nz, kx_grid_3D);
	Free_3D_mat(Nx, Ny, Nz, ky_grid_3D);
	Free_3D_mat(Nx, Ny, Nz, kz_grid_3D);	
			
	double complex *** in_film = AllocLargeArray_3D_Complex(Nx, Ny, Nz);
	for(int i = 0; i < Nx; i++){
	for(int j = 0; j < Ny; j++){
	for(int k = 0; k < Nz; k++){
		if(z_axis[k] > z_axis[sub_index] && z_axis[k] <= z_axis[film_index])
			in_film[i][j][k] = 1.0;
		else
			in_film[i][j][k] = 0.0;
	}}}
	
	double complex * in_film_RowMaj = malloc(sizeof(double complex) * Nx * Ny * Nz);
	Convert3D_2_RowMajor_Complex(Nx, Ny, Nz, in_film, in_film_RowMaj);
	
	Free_3D_mat_Complex(Nx, Ny, Nz, in_film);
	
	// --- Create Constants Structure --- //
	Constants -> Temperature = Temperature;
	Constants -> dt_factor = dt_factor;
	Constants -> ThermConst = ThermConst;
	Constants -> epsilon = epsilon;
	Constants -> BTO_pct = BTO_pct;
		
	Constants -> Nx = Nx; 
	Constants -> Ny = Ny;
	Constants -> Nz = Nz;

	Constants -> sub_index = sub_index;
	Constants -> film_index = film_index;
	Constants -> interface_index = interface_index;
	Constants -> Nz_film = Nz_film;
	
	Constants -> E_1_applied = E_1_applied;
	Constants -> E_2_applied = E_2_applied;
	Constants -> E_3_applied = E_3_applied;

	Constants -> Us_11 = Us_11;
	Constants -> Us_12 = Us_12;
	Constants -> Us_22 = Us_22;
	
	Constants -> LoadFromFile = LoadFromFile;
	Constants -> Pathname = Pathname;
	Constants -> Filename = malloc(strlen(Filename) + 1);
	strcpy(Constants -> Filename, Filename);
	
	Constants -> C11 = C11;
	Constants -> C12 = C12;
	Constants -> C44 = C44;
	
	Constants -> C = AllocLargeArray_4D(3, 3, 3, 3);
	// C is 0 indexed LOL!!
 	Constants -> C[0][0][0][0] = C11;
	Constants -> C[1][1][1][1] = C11;
	Constants -> C[2][2][2][2] = C11;
	
	Constants -> C[0][0][1][1] = C12;
	Constants -> C[0][0][2][2] = C12;
	Constants -> C[1][1][0][0] = C12;
	Constants -> C[1][1][2][2] = C12;
	Constants -> C[2][2][0][0] = C12;
	Constants -> C[2][2][1][1] = C12;
	
	Constants -> C[0][1][0][1] = C44;
	Constants -> C[1][0][1][0] = C44;
	Constants -> C[0][2][0][2] = C44;
	Constants -> C[2][0][2][0] = C44;
	Constants -> C[1][2][1][2] = C44;
	Constants -> C[2][1][2][1] = C44;
	Constants -> C[0][1][1][0] = C44;
	Constants -> C[1][0][0][1] = C44;
	Constants -> C[0][2][2][0] = C44;
	Constants -> C[2][0][0][2] = C44;
	Constants -> C[1][2][2][1] = C44;
	Constants -> C[2][1][1][2] = C44;

	Constants -> a1 = a1;
	Constants -> a11 = a11;
	Constants -> a12 = a12;
	Constants -> a111 = a111;
	Constants -> a112 = a112;
	Constants -> a123 = a123;
	Constants -> a1111 = a1111;
	Constants -> a1112 = a1112;
	Constants -> a1122 = a1122;
	Constants -> a1123 = a1123;
	
	Constants -> q11 = q11;
	Constants -> q12 = q12;
	Constants -> q44 = q44;
	
	Constants -> Q11 = Q11;
	Constants -> Q12 = Q12;
	Constants -> Q44 = Q44;
	
	Constants -> Lx = Lx;
	Constants -> Ly = Ly;
	Constants -> Lz = Lz;
	
	Constants -> G11 = G11;
	Constants -> G12 = G12;
	Constants -> H14 = H14;
	Constants -> H44 = H44;
	Constants -> G110 = G110;
	
	Constants -> dt = dt;
	Constants -> dt_factor = dt_factor;
	
	Constants -> k_electric_11 = k_electric_11;
	Constants -> k_electric_22 = k_electric_22;
	Constants -> k_electric_33 = k_electric_33;
	
	Constants -> saves = saves;
	
	Constants -> permittivity_0 = 8.85418782*1e-12;
	Constants -> boltzmann = 1.38064852 * 1e-23;
	
	Constants -> dx = dx;
	Constants -> dy = dy;
	Constants -> dz = dz;
	
	Constants -> x_axis = x_axis;
	Constants -> y_axis = y_axis;
	Constants -> z_axis = z_axis;

	Constants -> kx_axis = kx_axis;
	Constants -> ky_axis = ky_axis;
	Constants -> kz_axis = kz_axis;

	Constants -> h_film = z_axis[film_index];
	Constants -> h_sub = z_axis[0];
	Constants -> h_int = z_axis[interface_index];
	
	Constants -> kx_grid_3D_RowMaj = kx_grid_3D_RowMaj;
	Constants -> ky_grid_3D_RowMaj = ky_grid_3D_RowMaj;
	Constants -> kz_grid_3D_RowMaj = kz_grid_3D_RowMaj;

	Constants -> kx_grid_2D = kx_grid_2D;
	Constants -> ky_grid_2D = ky_grid_2D;

	Constants -> in_film_RowMaj = in_film_RowMaj;
	
	Constants -> NUM_THREADS = NUM_THREADS;
	Constants -> MULTITHREAD = MULTITHREAD;
	
	Constants -> kx_grid_2D_RowMaj = malloc(sizeof(double) * Nx * Ny);
	Constants -> ky_grid_2D_RowMaj = malloc(sizeof(double) * Nx * Ny);
	Convert2D_2_RowMajor(Nx, Ny, kx_grid_2D, Constants -> kx_grid_2D_RowMaj);
	Convert2D_2_RowMajor(Nx, Ny, ky_grid_2D, Constants -> ky_grid_2D_RowMaj);
	
	Constants -> total = Nx * Ny * Nz;
	Constants -> b11 = 0.5*C11*(Q11*Q11 + 2.0*Q12*Q12) + C12*Q12*(2.0*Q11 + Q12);
	Constants -> b12 = C11*Q12*(2.0*Q11 + Q12) + 
		C12*(Q11*Q11 + 3.0*Q12*Q12 + 2.0*Q11*Q12) + 2.0*C44*Q44*Q44;
	
	Constants -> MaxIter = MaxIter;
}

double coth( double x )
{
	return (1.0 / tanh(x));
}

// kspace axis setup
void Setup_k_axis(double L, int N, double * k)
{
	for(int i = 0; i < N; i++)
	{
		if(i <= N/2)
			k[i] = (2.0*M_PI/L) * i;
		else
			k[i] = (2.0*M_PI/L) * (i - N);
	}
}

void FreeSetup( struct ConstantsStruct * Constants )
{
	free(Constants -> x_axis);
	free(Constants -> y_axis);
	free(Constants -> z_axis);
	
	free(Constants -> kx_axis);
	free(Constants -> ky_axis);
	free(Constants -> kz_axis);
	
	free(Constants -> Filename);
	
	free(Constants -> kx_grid_3D_RowMaj);
	free(Constants -> ky_grid_3D_RowMaj);
	free(Constants -> kz_grid_3D_RowMaj);

	Free_2D_mat(Constants -> Nx, Constants -> Ny, Constants -> kx_grid_2D);
	Free_2D_mat(Constants -> Nx, Constants -> Ny, Constants -> ky_grid_2D);

	free(Constants -> kx_grid_2D_RowMaj);
	free(Constants -> ky_grid_2D_RowMaj);
	
	free(Constants -> in_film_RowMaj);
	
	Free_4D_mat(3, 3, 3, 3, Constants -> C);
	
	free(Constants -> saves);
	
}
