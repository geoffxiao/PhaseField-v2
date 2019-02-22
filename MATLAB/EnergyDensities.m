%%
% Constants Setup
Constants.dx = dx;
Constants.dy = dy;
Constants.dz = dz;

Constants.permittivity_0 = 8.85418782*1e-12;
permittivity_0 = Constants.permittivity_0;

ConstantsSetup;

% Set up real space axis
x_axis = linspace(0, Constants.Lx, Constants.Nx)';
y_axis = linspace(0, Constants.Ly, Constants.Ny)'; 
z_axis = linspace(0, Constants.Lz, Constants.Nz)';

% VERY IMPORTANT STEP!! Let us redefine zero, otherwise the exp() term will make our matrix solving go to shit...
z_axis = z_axis - z_axis((round(Constants.Nz/2)));

dx = x_axis(2) - x_axis(1); 
dy = y_axis(2) - y_axis(1); 
dz = z_axis(2) - z_axis(1);
[y_grid, x_grid, z_grid] = meshgrid(y_axis, x_axis, z_axis);
% kx_grid(x,y,z)

% Fourier space vectors -> Gradient energy terms and things that use air + film + sub
kx = 2*pi/Constants.Lx*[0:Constants.Nx/2 -Constants.Nx/2+1:-1]'; 
ky = 2*pi/Constants.Ly*[0:Constants.Ny/2 -Constants.Ny/2+1:-1]';
kz = 2*pi/Constants.Lz*[0:Constants.Nz/2 -Constants.Nz/2+1:-1]';
[ky_grid_3D,kx_grid_3D,kz_grid_3D] = meshgrid(ky,kx,kz);
[ky_grid_2D,kx_grid_2D] = meshgrid(ky,kx);

Constants.kx = kx;
Constants.ky = ky;
Constants.kz = kz;

Constants.x_axis = x_axis;
Constants.y_axis = y_axis;
Constants.z_axis = z_axis;    

Constants.dx = dx;
Constants.dy = dy;
Constants.dz = dz;
Constants.h_film = z_axis(Constants.film_index); % end of film
Constants.h_sub = z_axis(1); % where substrate ends...limit of elastic deformation allowed in substrate
Constants.h_int = z_axis(Constants.interface_index);

Constants.in_film = in_film;
ConstantsSetup;

%%
% Calculate everything
[strain_bc_mats_inv, strain_bc_mat_inv_korigin, ...
strain_sol_2Dk_1, strain_sol_2Dk_2, strain_sol_2Dk_3, ...
strain_sol_2Dk_4, strain_sol_2Dk_5, strain_sol_2Dk_6, ...
strain_sol_2Dk_d3_1, strain_sol_2Dk_d3_2, strain_sol_2Dk_d3_3, ...
strain_sol_2Dk_d3_4, strain_sol_2Dk_d3_5, strain_sol_2Dk_d3_6, ...
strain_bc_mats_L_inv_lu, strain_bc_mats_U_lu, strain_bc_mats_P_lu] = InfinitePlateSetup(Constants, kx_grid_2D, ky_grid_2D);

[electric_bc_mats_inv, electric_bc_mats_inv_korigin, ...
potential_2Dk_sol_1, potential_2Dk_sol_2, ...
potential_2Dk_d3_sol_1, potential_2Dk_d3_sol_2] = ElectrostaticSetup(Constants, kx_grid_2D, ky_grid_2D);       

[...
Green_11, Green_12, Green_13, ...
Green_21, Green_22, Green_23, ...
Green_31, Green_32, Green_33] = GreenTensorSetup(Constants, ...
                            kx_grid_3D, ky_grid_3D, kz_grid_3D);

[f1_0, f2_0, f3_0, u_1, u_2, u_3, Potential, ...
TotalStrain_11, TotalStrain_12, TotalStrain_13, ...
TotalStrain_22, TotalStrain_23, TotalStrain_33, ...
E_1_depol, E_2_depol, E_3_depol] = ...
                CalcEnergies(P1, P2, P3, Constants, ...
                        Green_11, Green_12, Green_13, ...
                        Green_21, Green_22, Green_23, ...
                        Green_31, Green_32, Green_33, ...
                        strain_bc_mats_inv, strain_bc_mat_inv_korigin, ...
                        electric_bc_mats_inv, electric_bc_mats_inv_korigin, ...
                        kx_grid_3D, ky_grid_3D, kz_grid_3D, z_axis, ...
                        strain_sol_2Dk_1, strain_sol_2Dk_2, strain_sol_2Dk_3, ...
                        strain_sol_2Dk_4, strain_sol_2Dk_5, strain_sol_2Dk_6, ...
                        strain_sol_2Dk_d3_1, strain_sol_2Dk_d3_2, strain_sol_2Dk_d3_3, ...
                        strain_sol_2Dk_d3_4, strain_sol_2Dk_d3_5, strain_sol_2Dk_d3_6, ...
                        potential_2Dk_sol_1, potential_2Dk_sol_2, ...
                        potential_2Dk_d3_sol_1, potential_2Dk_d3_sol_2, ...
                        strain_bc_mats_L_inv_lu, strain_bc_mats_U_lu, strain_bc_mats_P_lu);



%%
% Landau Energy
f_landau = a1 .* ( P1.^2 + P2.^2 + P3.^2 ) + a11 .* ( P1.^4 + P2.^4 + P3.^4 ) + ...
    a12 .* ( P1.^2 .* P2.^2 + P1.^2 .* P3.^2 + P2.^2 .* P3.^2 ) + ...
    a111 .* ( P1.^6 + P2.^6 + P3.^6 ) + a112 .* ( P1.^4 .* (P2.^2 + P3.^2) + P2.^4 .* (P3.^2 + P1.^2) + P3.^4 .* (P1.^2 + P2.^2) ) + ...
    a123 .* ( P1.^2 .* P2.^2 .* P3.^2 ) + ...
    a1111 .* ( P1.^8 + P2.^8 + P3.^8 ) + ...
    a1112 .* ( P1.^4 .* P2.^4 + P2.^4 .* P3.^4 + P1.^4 .* P3.^4 ) + ...
    a1123 .* ( P1.^4 .* P2.^2 .* P3.^2 + P2.^4 .* P3.^2 .* P1.^2 + P3.^4 .* P1.^2 .* P2.^2 );
f_l = f_landau;

%% 
% Elastic Energy
Eigenstrain_11 = Q11 * P1.^2 + Q12 .* (P2.^2 + P3.^2);
Eigenstrain_22 = Q11 * P2.^2 + Q12 .* (P1.^2 + P3.^2);
Eigenstrain_33 = Q11 * P3.^2 + Q12 .* (P1.^2 + P2.^2);
Eigenstrain_23 = Q44 * P2 .* P3;
Eigenstrain_13 = Q44 * P1 .* P3;
Eigenstrain_12 = Q44 * P1 .* P2; 

ElasticStrain_11 = TotalStrain_11 - Eigenstrain_11;
ElasticStrain_22 = TotalStrain_22 - Eigenstrain_22;
ElasticStrain_33 = TotalStrain_33 - Eigenstrain_33;
ElasticStrain_12 = TotalStrain_12 - Eigenstrain_12;
ElasticStrain_23 = TotalStrain_23 - Eigenstrain_23;
ElasticStrain_13 = TotalStrain_13 - Eigenstrain_13;

% Y L Li Elastic formulation
f_elastic = (C11/2) .* ( ElasticStrain_11.^2 + ElasticStrain_22.^2 + ElasticStrain_33.^2 ) + ...
    (C12) .* ( ElasticStrain_11 .* ElasticStrain_22 + ElasticStrain_22 .* ElasticStrain_33 + ElasticStrain_11 .* ElasticStrain_33 ) + ...
    (2*C44) .* ( ElasticStrain_12.^2 + ElasticStrain_23.^2 + ElasticStrain_13.^2 );

%%
% Het strain
e11_het = TotalStrain_11 - TotalStrain_homo_11(P1,P2,P3,Constants);
e22_het = TotalStrain_22 - TotalStrain_homo_22(P1,P2,P3,Constants);
e33_het = TotalStrain_33 - TotalStrain_homo_33(P1,P2,P3,Constants);
e12_het = TotalStrain_12 - TotalStrain_homo_12(P1,P2,P3,Constants);
e23_het = TotalStrain_23 - TotalStrain_homo_23(P1,P2,P3,Constants);
e13_het = TotalStrain_13 - TotalStrain_homo_13(P1,P2,P3,Constants);

e11_hom = TotalStrain_homo_11(P1,P2,P3,Constants);
e22_hom = TotalStrain_homo_22(P1,P2,P3,Constants);
e33_hom = TotalStrain_homo_33(P1,P2,P3,Constants);
e12_hom = TotalStrain_homo_12(P1,P2,P3,Constants);
e23_hom = TotalStrain_homo_23(P1,P2,P3,Constants);
e13_hom = TotalStrain_homo_13(P1,P2,P3,Constants);

%%
% Split Y L Li elastic into strain energy and electrostriction

% elastic + electrostriction
f_el = (C11/2) .* ( TotalStrain_11.^2 + TotalStrain_22.^2 + TotalStrain_33.^2 ) + ...
    (C12) .* ( TotalStrain_11 .* TotalStrain_22 + TotalStrain_22 .* TotalStrain_33 + TotalStrain_11 .* TotalStrain_33 ) + ...
    (2*C44) .* ( TotalStrain_12.^2 + TotalStrain_23.^2 + TotalStrain_13.^2 );
f_es = ...
    -( q11 .* TotalStrain_11 + q12 .* TotalStrain_22 + q12 .* TotalStrain_33 ) .* P1.^2 + ...
    -( q11 .* TotalStrain_22 + q12 .* TotalStrain_11 + q12 .* TotalStrain_33 ) .* P2.^2 + ...
    -( q11 .* TotalStrain_33 + q12 .* TotalStrain_11 + q12 .* TotalStrain_22 ) .* P3.^2 + ...
    -(2*q44) .* ( TotalStrain_12 .* P1 .* P2 + TotalStrain_23 .* P2 .* P3 + TotalStrain_13 .* P1 .* P3 );
b11 = 0.5*C11*(Q11^2 + 2*Q12^2) + C12*Q12*(2*Q11 + Q12);
b12 = C11*Q12*(2*Q11 + Q12) + C12*(Q11^2 + 3*Q12^2 + 2*Q11*Q12) + 2*C44*Q44^2;
f_prime_e = b11 .* ( P1.^4 + P2.^4 + P3.^4 ) + b12 .* ( P1.^2 .* P2.^2 + P1.^2 .* P3.^2 + P2.^2 .* P3.^2 );

%% Split into het and hom energies
% het
f_el_het = (C11/2) .* ( e11_het.^2 + e22_het.^2 + e33_het.^2 ) + ...
    (C12) .* ( e11_het .* e22_het + e22_het .* e33_het + e11_het .* e33_het ) + ...
    (2*C44) .* ( e12_het.^2 + e23_het.^2 + e13_het.^2 );
f_es_het = ...
    -( q11 .* e11_het + q12 .* e22_het + q12 .* e33_het ) .* P1.^2 + ...
    -( q11 .* e22_het + q12 .* e11_het + q12 .* e33_het ) .* P2.^2 + ...
    -( q11 .* e33_het + q12 .* e11_het + q12 .* e22_het ) .* P3.^2 + ...
    -(2*q44) .* ( e12_het .* P1 .* P2 + e23_het .* P2 .* P3 + e13_het .* P1 .* P3 );

% hom
f_el_hom = (C11/2) .* ( e11_hom.^2 + e22_hom.^2 + e33_hom.^2 ) + ...
    (C12) .* ( e11_hom .* e22_hom + e22_hom .* e33_hom + e11_hom .* e33_hom ) + ...
    (2*C44) .* ( e12_hom.^2 + e23_hom.^2 + e13_hom.^2 );
f_es_hom = ...
    -( q11 .* e11_hom + q12 .* e22_hom + q12 .* e33_hom ) .* P1.^2 + ...
    -( q11 .* e22_hom + q12 .* e11_hom + q12 .* e33_hom ) .* P2.^2 + ...
    -( q11 .* e33_hom + q12 .* e11_hom + q12 .* e22_hom ) .* P3.^2 + ...
    -(2*q44) .* ( e12_hom .* P1 .* P2 + e23_hom .* P2 .* P3 + e13_hom .* P1 .* P3 );

%%
% Gradient energy
P1_3Dk = fftn(P1); P2_3Dk = fftn(P2); P3_3Dk = fftn(P3);

P1_1 = real(ifftn(1i * P1_3Dk .* kx_grid_3D));
P1_2 = real(ifftn(1i * P1_3Dk .* ky_grid_3D));
P1_3 = real(ifftn(1i * P1_3Dk .* kz_grid_3D));

P2_1 = real(ifftn(1i * P2_3Dk .* kx_grid_3D));
P2_2 = real(ifftn(1i * P2_3Dk .* ky_grid_3D));
P2_3 = real(ifftn(1i * P2_3Dk .* kz_grid_3D));

P3_1 = real(ifftn(1i * P3_3Dk .* kx_grid_3D));
P3_2 = real(ifftn(1i * P3_3Dk .* ky_grid_3D));
P3_3 = real(ifftn(1i * P3_3Dk .* kz_grid_3D));

f_grad = (G11/2) * ( P1_1.^2 + P2_2.^2 + P3_3.^2 ) + ...
    G12 .* ( P1_1 .* P2_2 + P2_2 .* P3_3 + P3_3 .* P1_1 ) + ...
    (H44/2) .* ( P1_2.^2 + P2_1.^2 + P2_3.^2 + P3_2.^2 + P1_3.^2 + P3_1.^2 ) + ...
    (H14 - G12) .* ( P1_2 .* P2_1 + P1_3 .* P3_1 + P2_3 .* P3_2 );

%%
% Electric energy
permittivity_0 = Constants.permittivity_0;
f_elec = - (1/2) * ( E_1_depol .* P1 + E_2_depol .* P2 + E_3_depol .* P3 ) + ...
    (-1/2) * permittivity_0 * (k_electric_11 * E_1_depol.^2 + k_electric_22 * E_2_depol.^2 + k_electric_33 * E_3_depol.^2);

f_elec_ext = -( E_1_applied .* P1 + E_2_applied .* P2 + E_3_applied .* P3 );

%%
% total energy
f = f_l + f_el + f_es + f_elec + f_grad + f_prime_e;


%% 
% total energy
f_landau = f_landau .* in_film;
f_elastic = f_elastic .* in_film;
f_grad = f_grad .* in_film;
f_elec = f_elec .* in_film;
f_elec_ext = f_elec_ext .* in_film;
f_tot=f_landau+f_elastic+f_grad+f_elec;

% F_landau = trapz(trapz(trapz(f_landau)));
% F_elastic = trapz(trapz(trapz(f_elastic)));
% F_grad = trapz(trapz(trapz(f_grad)));
% F_elec = trapz(trapz(trapz(f_elec)));
% F_elec_ext = trapz(trapz(trapz(f_elec_ext)));

% F_landau = sum(sum(sum(f_landau)));
% F_elastic = sum(sum(sum(f_elastic)));
% F_grad = sum(sum(sum(f_grad)));
% F_elec = sum(sum(sum(f_elec)));
% F_elec_ext = sum(sum(sum(f_elec_ext)));

F_landau = trapz(trapz(trapz(f_landau,3),2),1);
F_elastic = trapz(trapz(trapz(f_elastic,3),2),1);
F_grad = trapz(trapz(trapz(f_grad,3),2),1);
F_elec = trapz(trapz(trapz(f_elec,3),2),1);
F_elec_ext = trapz(trapz(trapz(f_elec_ext,3),2),1);


[F_landau, F_elastic, F_grad, F_elec]

F_tot = sum([F_landau, F_elastic,F_grad,F_elec]);