Temperature = 25;
dt_factor = 0.03;
BTO_pct = 0.6; 
l_0 = 1e-9; 
l_x = 1e-9; l_y = 1e-9; l_z = 1e-9;
Nx = 48; Ny = 48; Nz = 32; % Grid Points
GradientEnergyFactor = 1;

ThermConst = 0;
epsilon = 1e-3; % convergence criterion
saves = [0 : 100 : 5000]; % save after this many iterations

Substrates;
Substrate_lattice_11 = SSO_lattice_11;
Substrate_lattice_22 = SSO_lattice_22;

sub_index = 8; % where substrate starts, >= 1, at and below this index P = 0, substrate thickness
film_index = Nz - 2; % where film ends, <= Nz, above this index, P = 0

LOAD = 1; % 1 = Load initial conditions from file called init.mat, 0 = random P initial conditions


%% Grid Size
% sub_index = 0 to get no substrate part
interface_index = sub_index + 1;
Nz_film = film_index - sub_index; % film lies in the area between sub_index + 1 and film_index, film thickness
Nz_film_sub = film_index;

%%
% Run inputs
rng(1)
E_1_applied = 0; E_2_applied = 0; E_3_applied = 0; % in V/m, 1e5 V/m = 1 kV/cm

% http://www.mmm.psu.edu/JWang2010_JAP_BTO1stPrinciple.pdf
BTO_lattice = 4.006;
% https://deepblue.lib.umich.edu/bitstream/handle/2027.42/62658/nature02773.pdf?sequence=1&isAllowed=y
STO_lattice = 3.905;
STO_pct = 1 - BTO_pct;
BST_lattice = BTO_pct * BTO_lattice + STO_pct * STO_lattice;

% http://www.iue.tuwien.ac.at/phd/dhar/node12.html
Us_11 = ((Substrate_lattice_11 - BST_lattice) / BST_lattice); % unitless misfit strain
Us_22 = ((Substrate_lattice_22 - BST_lattice) / BST_lattice); % anisotropic misfit strain
Us_12 = 0;

PATH = './';
STRING = sprintf('BTO_%g-%g',BTO_pct*1e2,STO_pct*1e2); % Material Name

VPA_ELECTRIC_ON = 0; % numerical errors when doing electric energy, so we need to use vpa
VPA_ELASTIC_ON = 0;

%% ---- Nothing needed to modify below ---- %%
%% Energies, 0 = off, 1 = on
ELASTIC = 1;
HET_ELASTIC_RELAX = 1 && ELASTIC; % automatically 0 if no elastic energy 
% no heterogenous elastic relaxation if no elastic energy

ELECTRIC = 1;
%% BTO =  http://www.ems.psu.edu/~chen/publications/YL2005APL.pdf
C11_BTO = 1.78 * 1e11;
C12_BTO = 0.964 * 1e11;
C44_BTO = 1.22 * 1e11;

Q11_BTO = 0.10; 
Q12_BTO = -0.034;
Q44_BTO = 0.029;

a1_BTO = @(T) 4.124 * ( T - 115 ) * 1e5;
a11_BTO = -2.097 * 1e8;
a12_BTO = 7.974 * 1e8;
a111_BTO = 1.294 * 1e9;
a112_BTO = -1.950 * 1e9;
a123_BTO = -2.5 * 1e9;
a1111_BTO = 3.863 * 1e10;
a1112_BTO = 2.529 * 1e10;
a1122_BTO = 1.637 * 1e10;
a1123_BTO = 1.367 * 1e10;

%% STO, Table VII, Table I last column
% http://www.wolframalpha.com/input/?i=1+dyn%2Fcm%5E2
% 1 dyn/cm^2 = 0.1 Pa
C11_STO = 3.36 * 1e11;
C12_STO = 1.07 * 1e11;
C44_STO = 1.27 * 1e11;

% STO Yu Luan Li PRB 184112, Table II
% 1 cm^4 / esu^2 = 8.988 * 1e10 m^4 / C^2
% http://www.wolframalpha.com/input/?i=(1+cm)%5E4%2F(3.3356e-10+coulomb)%5E2
Q11_STO = 4.57 * 1e-2; 
Q12_STO = -1.348 * 1e-2;
Q44_STO = 0.957 * 1e-2;

% STO Yu Luan Li PRB 184112, Table III a1 (1st, from ref 6, 16, 22), Table
% VII a11, a12, a1
% Shirokov, Yuzyu, PRB 144118, Table I STO in parantheses

% 1 cm^2 * dyn / esu^2 = 8.9876 * 1e9 J * m^2 / C
a1_STO = @(T) 4.05 * 1e7 * ( coth( 54 / (T+273) ) - coth(54/30) );

% 1 cm^6 * dyn / esu^4 = 8.0776 * 1e20 J * m^5 / C^4
a11_STO = 17 * 1e8;
a12_STO = 13.7 * 1e8;

a111_STO = 0;
a112_STO = 0;
a123_STO = 0;
a1111_STO = 0;
a1112_STO = 0;
a1122_STO = 0;
a1123_STO = 0;

%% Elastic Tensor
C11 = C11_BTO * BTO_pct + C11_STO * STO_pct;
C12 = C12_BTO * BTO_pct + C12_STO * STO_pct;
C44 = C44_BTO * BTO_pct + C44_STO * STO_pct;

C = zeros(3,3,3,3);
C(1,1,1,1) = C11; C(2,2,2,2) = C11; C(3,3,3,3) = C11;
C(1,1,2,2) = C12; C(1,1,3,3) = C12; C(2,2,1,1) = C12;
C(2,2,3,3) = C12; C(3,3,1,1) = C12; C(3,3,2,2) = C12;
C(1,2,1,2) = C44; C(2,1,2,1) = C44; C(1,3,1,3) = C44;
C(3,1,3,1) = C44; C(2,3,2,3) = C44; C(3,2,3,2) = C44;
C(1,2,2,1) = C44; C(2,1,1,2) = C44; C(1,3,3,1) = C44;
C(3,1,1,3) = C44; C(2,3,3,2) = C44; C(3,2,2,3) = C44;

%% Electrostriction
q11_BTO = C11_BTO * Q11_BTO + 2 * C12_BTO * Q12_BTO;
q12_BTO = C11_BTO * Q12_BTO + C12_BTO * (Q11_BTO + Q12_BTO);
q44_BTO = 2 * C44_BTO * Q44_BTO;

q11_STO = C11_STO * Q11_STO + 2 * C12_STO * Q12_STO;
q12_STO = C11_STO * Q12_STO + C12_STO * (Q11_STO + Q12_STO);
q44_STO = 2 * C44_STO * Q44_STO;

q11 = q11_BTO * BTO_pct + q11_STO * STO_pct;
q12 = q12_BTO * BTO_pct + q12_STO * STO_pct;
q44 = q44_BTO * BTO_pct + q44_STO * STO_pct;

q11_h = q11 + 2 * q12;
q22_h = q11 - q12;
C11_h = C11 + 2 * C12;
C22_h = C11 - C12;

Q11 = (1/3) * ( (q11_h/C11_h) + (2 * q22_h/C22_h) );
Q12 = (1/3) * ( (q11_h/C11_h) - (q22_h/C22_h) );
Q44 = q44 / (2 * C44);

%% LGD Constants, stress free BTO
a1_T = @(T) a1_BTO(T) * BTO_pct + a1_STO(T) * STO_pct;
a1 = a1_T(Temperature);
a11 = a11_BTO * BTO_pct + a11_STO * STO_pct;
a12 = a12_BTO * BTO_pct + a12_STO * STO_pct;
a111 = a111_BTO * BTO_pct + a111_STO * STO_pct;
a112 = a112_BTO * BTO_pct + a112_STO * STO_pct;
a123 = a123_BTO * BTO_pct + a123_STO * STO_pct;
a1111 = a1111_BTO * BTO_pct + a1111_STO * STO_pct;
a1112 = a1112_BTO * BTO_pct + a1112_STO * STO_pct;
a1122 = a1122_BTO * BTO_pct + a1122_STO * STO_pct;
a1123 = a1123_BTO * BTO_pct + a1123_STO * STO_pct;

%% Simulation Size
Lx = l_x*(Nx-1); Ly = l_y*(Ny-1); Lz = l_z*(Nz-1);

%% Gradient Free Energy
% Let us just make all BST comps have the same gradient energy...
a0 = abs(a1_BTO(25));
G110 = l_0^2 * a0; 
G11 = GradientEnergyFactor * G110;

G12 = 0; % anisotropic part

% G44 = G44'
H14 = 0; % H14 = G12 + G44 - G44'
H44 = G11; % H44 = G44 + G44'

%% Run time setup
t_0 = 1/a0;
dt = dt_factor/(a0);
RUN_TIME = 1; % time, not iterations
MAX_ITERATIONS = RUN_TIME / dt;

%% Electrical
k_electric_11 = 100; k_electric_22 = k_electric_11; k_electric_33 = k_electric_11;

%% Constants struct
Constants.E_1_applied = E_1_applied;
Constants.E_2_applied = E_2_applied;
Constants.E_3_applied = E_3_applied;
Constants.Us_11 = Us_11;
Constants.Us_22 = Us_22;
Constants.Us_12 = Us_12;
Constants.Temperature = Temperature;
Constants.Nx = Nx;
Constants.Ny = Ny;
Constants.Nz = Nz;
Constants.sub_index = sub_index;
Constants.interface_index = interface_index;
Constants.film_index = film_index;
Constants.ELASTIC = ELASTIC;
Constants.HET_ELASTIC_RELAX = HET_ELASTIC_RELAX;
Constants.ELECTRIC = ELECTRIC;
Constants.ThermConst = ThermConst;
Constants.C11 = C11;
Constants.C12 = C12;
Constants.C44 = C44;
Constants.C = C;
Constants.q11 = q11;
Constants.q12 = q12;
Constants.q44 = q44;
Constants.Q11 = Q11;
Constants.Q12 = Q12;
Constants.Q44 = Q44;
Constants.a1_T = a1_T;
Constants.a1 = a1;
Constants.a11 = a11;
Constants.a12 = a12;
Constants.a111 = a111;
Constants.a112 = a112;
Constants.a123 = a123;
Constants.a1111 = a1111;
Constants.a1112 = a1112;
Constants.a1122 = a1122;
Constants.a1123 = a1123;
Constants.Lx = Lx;
Constants.Ly = Ly;
Constants.Lz = Lz;
Constants.G11 = G11;
Constants.G12 = G12;
Constants.H14 = H14;
Constants.H44 = H44;
Constants.dt = dt;
Constants.k_electric_11 = k_electric_11;
Constants.k_electric_22 = k_electric_22;
Constants.k_electric_33 = k_electric_33;
Constants.BTO_pct = BTO_pct;
Constants.saves = saves;
Constants.permittivity_0 = 8.85418782*1e-12;
Constants.boltzmann = 1.38064852 * 1e-23;    

%%
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

%
if( sub_index > 0 )
    in_film = (z_grid > z_axis(sub_index)) & (z_grid <= z_axis(film_index));
    not_in_film = ~in_film;
    Constants.not_in_film = +not_in_film;
    Constants.in_film = +in_film;
else
    Constants.in_film = 1;
    Constants.not_in_film = 0;
end

% z axis
Constants.h_film = z_axis(Constants.film_index); % end of film
Constants.h_sub = z_axis(1); % where substrate ends...limit of elastic deformation allowed in substrate
Constants.h_int = z_axis(Constants.interface_index);
