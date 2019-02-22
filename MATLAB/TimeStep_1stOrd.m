function [P1, P2, P3] = TimeStep_1stOrd(P1_prev, P2_prev, P3_prev, f1_prev, f2_prev, f3_prev, ...
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
    
    P1_prev_3Dk = fftn(P1_prev);
    P2_prev_3Dk = fftn(P2_prev);
    P3_prev_3Dk = fftn(P3_prev);
    
    f1_prev_3Dk = fftn(f1_prev);
    f2_prev_3Dk = fftn(f2_prev);
    f3_prev_3Dk = fftn(f3_prev);
    
    % Time step
    P1_A_3Dk = (( P1_prev_3Dk + dt .* -f1_prev_3Dk ) ./ ...
                ( 1 + dt .* ( H44 .* (ky_grid_3D.^2 + kz_grid_3D.^2 ) + (G11 .* kx_grid_3D.^2) ) ));
    P2_A_3Dk = (( P2_prev_3Dk + dt .* -f2_prev_3Dk ) ./ ...
                ( 1 + dt .* ( H44 .* (kx_grid_3D.^2 + kz_grid_3D.^2 ) + (G11 .* ky_grid_3D.^2) ) ));
    P3_A_3Dk = (( P3_prev_3Dk + dt .* -f3_prev_3Dk ) ./ ...
                ( 1 + dt .* ( H44 .* (kx_grid_3D.^2 + ky_grid_3D.^2 ) + (G11 .* kz_grid_3D.^2) ) ));

    P1_A = real(ifftn(P1_A_3Dk));
    P2_A = real(ifftn(P2_A_3Dk));
    P3_A = real(ifftn(P3_A_3Dk));

    P1_A_2Dk = fft_2d_slices(P1_A);
    P2_A_2Dk = fft_2d_slices(P2_A);
    P3_A_2Dk = fft_2d_slices(P3_A);

    Time_bc_1 = cat(3, -P1_A_2Dk(:,:,interface_index), -P1_A_2Dk(:,:,film_index));
    Time_bc_2 = cat(3, -P2_A_2Dk(:,:,interface_index), -P2_A_2Dk(:,:,film_index));
    Time_bc_3 = cat(3, -P3_A_2Dk(:,:,interface_index), -P3_A_2Dk(:,:,film_index));

    Time_C_1 = zeros(Nx, Ny, 2);
    Time_C_2 = zeros(Nx, Ny, 2);
    Time_C_3 = zeros(Nx, Ny, 2);

    for k1_ind = 1 : Nx
    for k2_ind = 1 : Ny

        Time_C_1(k1_ind,k2_ind,:) = squeeze(TimeCoeffs_inv_mat_1(k1_ind,k2_ind,:,:)) * ...
                   squeeze(Time_bc_1(k1_ind,k2_ind,:));

        Time_C_2(k1_ind,k2_ind,:) = squeeze(TimeCoeffs_inv_mat_2(k1_ind,k2_ind,:,:)) * ...
                   squeeze(Time_bc_2(k1_ind,k2_ind,:));

        Time_C_3(k1_ind,k2_ind,:) = squeeze(TimeCoeffs_inv_mat_3(k1_ind,k2_ind,:,:)) * ...
                   squeeze(Time_bc_3(k1_ind,k2_ind,:));

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

    P1 = P1_A + P1_B;
    P2 = P2_A + P2_B;
    P3 = P3_A + P3_B;

    P1 = P1 .* in_film;
    P2 = P2 .* in_film;
    P3 = P3 .* in_film;

end