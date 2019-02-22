function [electric_bc_mats_inv, electric_bc_mats_inv_korigin] = ElectrostaticSetup_vpa(Constants, kx_grid_2D, ky_grid_2D)
    
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    k_electric_11 = Constants.k_electric_11;
    k_electric_22 = Constants.k_electric_22;
    k_electric_33 = Constants.k_electric_33;
    h_int = Constants.h_int;
    h_film = Constants.h_film;
    
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
     
end