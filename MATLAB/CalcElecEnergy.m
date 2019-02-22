function [f1_elec, f2_elec, f3_elec, Potential, ...
            E_1_depol, E_2_depol, E_3_depol] = CalcElecEnergy(P1, P2, P3, Constants, ...
                                                    electric_bc_mats_inv, electric_bc_mats_inv_korigin, ...
                                                    kx_grid_3D, ky_grid_3D, kz_grid_3D, z_axis, ...
                                                    potential_2Dk_sol_1, potential_2Dk_sol_2, ...
                                                    potential_2Dk_d3_sol_1, potential_2Dk_d3_sol_2)

    k_electric_11 = Constants.k_electric_11;
    k_electric_22 = Constants.k_electric_22;
    k_electric_33 = Constants.k_electric_33;
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    Nz = Constants.Nz;
    interface_index = Constants.interface_index;
    film_index = Constants.film_index;
    permittivity_0 = Constants.permittivity_0;
    ELECTRIC = Constants.ELECTRIC;

    k_mag_3D = sqrt(kx_grid_3D.^2 + ky_grid_3D.^2 + kz_grid_3D.^2);
    kx_grid_2D = kx_grid_3D(:,:,1);
    ky_grid_2D = ky_grid_3D(:,:,1);  
    
    if( ELECTRIC )

         %% Solve particular solution within the entire sample...
        P1_3Dk = fftn(P1); P2_3Dk = fftn(P2); P3_3Dk = fftn(P3);

        Potential_A_k = -1i .* (kx_grid_3D .* P1_3Dk + ky_grid_3D .* P2_3Dk + kz_grid_3D .* P3_3Dk) ./ ...
            ( permittivity_0 .* (k_electric_11 .* kx_grid_3D.^2 + k_electric_22 .* ky_grid_3D.^2 + k_electric_33 .* kz_grid_3D.^2) );
        Potential_A_k(k_mag_3D == 0) = 0;
        Potential_A = real(ifftn(Potential_A_k));

        E_A_1_k = -1i .* kx_grid_3D .* Potential_A_k;
        E_A_2_k = -1i .* ky_grid_3D .* Potential_A_k;
        E_A_3_k = -1i .* kz_grid_3D .* Potential_A_k;

        E_A_1 = real(ifftn(E_A_1_k));
        E_A_2 = real(ifftn(E_A_2_k));
        E_A_3 = real(ifftn(E_A_3_k));

        %% Homogenous solution boundary condition application, 2D FT solution
        C_mat = zeros(Nx,Ny,2);
        potential_bc_interface = squeeze(Potential_A(:,:,interface_index));
        potential_bc_film = squeeze(Potential_A(:,:,film_index));
        potential_bc_interface_2Dk = -fft2(potential_bc_interface);
        potential_bc_film_2Dk = -fft2(potential_bc_film);
        potential_bc_given_mat = cat(3,potential_bc_interface_2Dk,potential_bc_film_2Dk);

        Potential_B_2Dk = zeros(Nx,Ny,Nz);
        Potential_B_2Dk_d3 = zeros(Nx,Ny,Nz);

        temp1 = zeros(2,2);
        temp2 = zeros(2,1);
        
        for k1 = 1 : Nx
        for k2 = 1 : Ny
            temp1(:) = electric_bc_mats_inv(k1,k2,:,:);
            temp2(:) = potential_bc_given_mat(k1,k2,:);
%            C_mat(k1,k2,:) = squeeze(electric_bc_mats_inv(k1,k2,:,:)) * squeeze(potential_bc_given_mat(k1,k2,:));
            C_mat(k1,k2,:) = temp1 * temp2;
        end
        end

        C1 = squeeze(C_mat(:,:,1)); C2 = squeeze(C_mat(:,:,2));
        C1 = repmat(C1,1,1,Nz);
        C2 = repmat(C2,1,1,Nz);
        p = sqrt((k_electric_11*kx_grid_2D.^2 + k_electric_22*ky_grid_2D.^2)/k_electric_33);
        

        Potential_B_2Dk = C1 .* potential_2Dk_sol_1 + C2 .* potential_2Dk_sol_2;
        Potential_B_2Dk_d3 = C1 .* potential_2Dk_d3_sol_1 + C2 .* potential_2Dk_d3_sol_2;

        % Fourier space origin
        C_origin = electric_bc_mats_inv_korigin * squeeze(potential_bc_given_mat(1,1,:));
        C1 = C_origin(1); C2 = C_origin(2);
        Potential_B_2Dk(1,1,:) = C1 * z_axis + C2;
        Potential_B_2Dk_d3(1,1,:) = C1;

        %% Total Potential
        Potential_B = ifft_2d_slices(Potential_B_2Dk);
        Potential = Potential_A + Potential_B;

        % E field
        E_B_2Dk_1 = -1i .* kx_grid_3D .* Potential_B_2Dk;
        E_B_2Dk_2 = -1i .* ky_grid_3D .* Potential_B_2Dk;
        E_B_2Dk_3 = -Potential_B_2Dk_d3;

        E_B_1 = ifft_2d_slices(E_B_2Dk_1);
        E_B_2 = ifft_2d_slices(E_B_2Dk_2);
        E_B_3 = ifft_2d_slices(E_B_2Dk_3);

        E_1_depol = E_A_1 + E_B_1;
        E_2_depol = E_A_2 + E_B_2;
        E_3_depol = E_A_3 + E_B_3;

        %% Electrical energy define
        f1_elec = -0.5 * E_1_depol; 
        f2_elec = -0.5 * E_2_depol; 
        f3_elec = -0.5 * E_3_depol; 
        
    else

        f1_elec = 0; f2_elec = 0; f3_elec = 0;
        Potential = 0;

    end
    
end