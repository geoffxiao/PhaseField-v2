function [P1_2, P2_2, P3_3] = TimeStep_2ndOrd(P1_0, P2_0, P3_0, ...
                            P1_1, P2_1, P3_1, ...
                            f1_0, f2_0, f3_0, ...
                            f1_1, f2_1, f3_1, ...
                            TimeCoeffs_inv_mat_1, TimeCoeffs_inv_mat_2, TimeCoeffs_inv_mat_3, ...
                            q_mat_1, q_mat_2, q_mat_3, ...
                            scaling_1, scaling_2, scaling_3, ...
                            Constants, kx_grid_3D, ky_grid_3D, kz_grid_3D, z_axis)
                        
    G11 = Constants.G11;
    H44 = Constants.H44;
    dt = Constants.dt;
    in_film = Constants.in_film;
    interface_index = Constants.interface_index;
    film_index = Constants.film_index;
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    Nz = Constants.Nz;
    
    P1_0_3Dk = fftn(P1_0);
    P2_0_3Dk = fftn(P2_0);
    P3_0_3Dk = fftn(P3_0);
    
    P1_1_3Dk = fftn(P1_1);
    P2_1_3Dk = fftn(P2_1);
    P3_1_3Dk = fftn(P3_1);    
    
    f1_0_3Dk = fftn(f1_0);
    f2_0_3Dk = fftn(f2_0);
    f3_0_3Dk = fftn(f3_0);
 
    f1_1_3Dk = fftn(f1_1);
    f2_1_3Dk = fftn(f2_1);
    f3_1_3Dk = fftn(f3_1);
    
    % Time step
    P1_2_A_3Dk = ((4 .* P1_1_3Dk - P1_0_3Dk - 2 .* dt .* (2 .* f1_1_3Dk - f1_0_3Dk)) ./ ...
                 (3 + 2 .* dt .* kx_grid_3D.^2 .* G11 + 2 .* dt .* H44 .* ky_grid_3D.^2 + 2 .* dt .* H44 .* kz_grid_3D.^2));
    
    P2_2_A_3Dk = ((4 .* P2_1_3Dk - P2_0_3Dk - 2 .* dt .* (2 .* f2_1_3Dk - f2_0_3Dk)) ./ ...
                 (3 + 2 .* dt .* ky_grid_3D.^2 .* G11 + 2 .* dt .* H44 .* kx_grid_3D.^2 + 2 .* dt .* H44 .* kz_grid_3D.^2));
        
    P3_2_A_3Dk = ((4 .* P3_1_3Dk - P3_0_3Dk - 2 .* dt .* (2 .* f3_1_3Dk - f3_0_3Dk)) ./ ...
                 (3 + 2 .* dt .* kz_grid_3D.^2 .* G11 + 2 .* dt .* H44 .* kx_grid_3D.^2 + 2 .* dt .* H44 .* ky_grid_3D.^2));
    
    P1_2_A = real(ifftn(P1_2_A_3Dk));
    P2_2_A = real(ifftn(P2_2_A_3Dk));
    P3_2_A = real(ifftn(P3_2_A_3Dk));

    P1_2_A_2Dk = fft_2d_slices(P1_2_A);
    P2_2_A_2Dk = fft_2d_slices(P2_2_A);
    P3_2_A_2Dk = fft_2d_slices(P3_2_A);

    Time_bc_1 = cat(3, -P1_2_A_2Dk(:,:,interface_index), -P1_2_A_2Dk(:,:,film_index));
    Time_bc_2 = cat(3, -P2_2_A_2Dk(:,:,interface_index), -P2_2_A_2Dk(:,:,film_index));
    Time_bc_3 = cat(3, -P3_2_A_2Dk(:,:,interface_index), -P3_2_A_2Dk(:,:,film_index));

    Time_C_1 = zeros(Nx, Ny, 2);
    Time_C_2 = zeros(Nx, Ny, 2);
    Time_C_3 = zeros(Nx, Ny, 2);

    temp1 = zeros(2,2);
    temp2 = zeros(2,1);
    for k1_ind = 1 : Nx
    for k2_ind = 1 : Ny

        temp1(:) = TimeCoeffs_inv_mat_1(k1_ind,k2_ind,:,:);
        temp2(:) = Time_bc_1(k1_ind,k2_ind,:);
        Time_C_1(k1_ind,k2_ind,:) = temp1 * temp2;

        temp1(:) = TimeCoeffs_inv_mat_2(k1_ind,k2_ind,:,:);
        temp2(:) = Time_bc_2(k1_ind,k2_ind,:);        
        Time_C_2(k1_ind,k2_ind,:) = temp1 * temp2;

        temp1(:) = TimeCoeffs_inv_mat_3(k1_ind,k2_ind,:,:);
        temp2(:) = Time_bc_3(k1_ind,k2_ind,:);
        Time_C_3(k1_ind,k2_ind,:) = temp1 * temp2;

    end
    end

    P1_B_2Dk = zeros(Nx,Ny,Nz);
    P2_B_2Dk = zeros(Nx,Ny,Nz);
    P3_B_2Dk = zeros(Nx,Ny,Nz);

    for z_loop = 1 : Nz
        A = Time_C_1(:,:,1); B = Time_C_1(:,:,2);
        q = q_mat_1;
        P1_B_2Dk(:,:,z_loop) = A .* exp(q * z_axis(z_loop) - scaling_1) + B .* exp(-q * z_axis(z_loop) - scaling_1);
    end

    for z_loop = 1 : Nz
        A = Time_C_2(:,:,1); B = Time_C_2(:,:,2);
        q = q_mat_2;
        P2_B_2Dk(:,:,z_loop) = A .* exp(q * z_axis(z_loop) - scaling_2) + B .* exp(-q * z_axis(z_loop) - scaling_2);
    end

    for z_loop = 1 : Nz
        A = Time_C_3(:,:,1); B = Time_C_3(:,:,2);
        q = q_mat_3;
        P3_B_2Dk(:,:,z_loop) = A .* exp(q * z_axis(z_loop) - scaling_3) + B .* exp(-q * z_axis(z_loop) - scaling_3);
    end

    P1_B = ifft_2d_slices( P1_B_2Dk );
    P2_B = ifft_2d_slices( P2_B_2Dk );
    P3_B = ifft_2d_slices( P3_B_2Dk );

    P1_2 = P1_2_A + P1_B;
    P2_2 = P2_2_A + P2_B;
    P3_3 = P3_2_A + P3_B;

    P1_2 = P1_2 .* in_film;
    P2_2 = P2_2 .* in_film;
    P3_3 = P3_3 .* in_film;

end