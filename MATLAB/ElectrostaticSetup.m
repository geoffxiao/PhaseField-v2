function [electric_bc_mats_inv, electric_bc_mats_inv_korigin, ...
    potential_2Dk_sol_1, potential_2Dk_sol_2, ...
    potential_2Dk_d3_sol_1, potential_2Dk_d3_sol_2] = ElectrostaticSetup(Constants, kx_grid_2D, ky_grid_2D)
    
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    k_electric_11 = Constants.k_electric_11;
    k_electric_22 = Constants.k_electric_22;
    k_electric_33 = Constants.k_electric_33;
    h_int = Constants.h_int;
    h_film = Constants.h_film;
    Nz = Constants.Nz;
    z_axis = Constants.z_axis;

    electric_bc_mats_inv = zeros(Nx,Ny,2,2);

    for k1_ind = 1 : Nx
    for k2_ind = 1 : Ny
        if(k1_ind ~= 1 || k2_ind ~= 1)
            eta_1 = kx_grid_2D(k1_ind,k2_ind);
            eta_2 = ky_grid_2D(k1_ind,k2_ind);
            p = sqrt((k_electric_11*eta_1^2 + k_electric_22*eta_2^2)/k_electric_33);

            potential_sol_mat = [ exp(h_int*p), exp(-h_int*p); ...
                                  exp(h_film*p), exp(-h_film*p)];

            electric_bc_mats_inv(k1_ind,k2_ind,:,:) = ...
                1 / (potential_sol_mat(1,1) * potential_sol_mat(2,2) - ...
                     potential_sol_mat(1,2) * potential_sol_mat(2,1)) * ...
                [ potential_sol_mat(2,2) -potential_sol_mat(1,2); ...
                 -potential_sol_mat(2,1) potential_sol_mat(1,1)];
        end            
    end
    end

    electric_bc_mats_inv_korigin = [h_int,  1;...
                                    h_film, 1];
    electric_bc_mats_inv_korigin = ...
        1 / (electric_bc_mats_inv_korigin(1,1) * electric_bc_mats_inv_korigin(2,2) - ...
             electric_bc_mats_inv_korigin(1,2) * electric_bc_mats_inv_korigin(2,1)) * ...
        [ electric_bc_mats_inv_korigin(2,2) -electric_bc_mats_inv_korigin(1,2); ...
         -electric_bc_mats_inv_korigin(2,1) electric_bc_mats_inv_korigin(1,1)];

    potential_2Dk_sol_1 = zeros(Nx,Ny,Nz);
    potential_2Dk_sol_2 = zeros(Nx,Ny,Nz);
    potential_2Dk_d3_sol_1 = zeros(Nx,Ny,Nz);
    potential_2Dk_d3_sol_2 = zeros(Nx,Ny,Nz);

    p = sqrt((k_electric_11*kx_grid_2D.^2 + k_electric_22*ky_grid_2D.^2)/k_electric_33);
    for z_loop = 1 : Nz
        potential_2Dk_sol_1(:,:,z_loop) = exp(z_axis(z_loop) .* p);
        potential_2Dk_sol_2(:,:,z_loop) = exp(-z_axis(z_loop) .* p);
        potential_2Dk_d3_sol_1(:,:,z_loop) = p .* exp(z_axis(z_loop) .* p);
        potential_2Dk_d3_sol_2(:,:,z_loop) = -p .* exp(-z_axis(z_loop) .* p);
    end     
end