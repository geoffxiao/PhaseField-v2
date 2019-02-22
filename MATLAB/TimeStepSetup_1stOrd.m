function [q_mat_1, q_mat_2, q_mat_3, ...
    scaling_1, scaling_2, scaling_3, ...
    TimeCoeffs_inv_mat_1, TimeCoeffs_inv_mat_2, TimeCoeffs_inv_mat_3] = ...
    TimeStepSetup_1stOrd(Constants, kx_grid_2D, ky_grid_2D, z_axis)

    dt = Constants.dt;
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    interface_index = Constants.interface_index;
    h_int = Constants.h_int;
    film_index = Constants.film_index;
    h_film = Constants.h_film;
    G11 = Constants.G11;
    H44 = Constants.H44;

    TimeCoeffs_inv_mat_1 = zeros(Nx, Ny, 2, 2);
    TimeCoeffs_inv_mat_2 = zeros(Nx, Ny, 2, 2);
    TimeCoeffs_inv_mat_3 = zeros(Nx, Ny, 2, 2);

    q_mat_1 = zeros(Nx, Ny);
    q_mat_2 = zeros(Nx, Ny);
    q_mat_3 = zeros(Nx, Ny);

    scaling_1 = zeros(Nx,Ny);
    scaling_2 = zeros(Nx,Ny);
    scaling_3 = zeros(Nx,Ny);

    % BC at h_int
    % BC at h_film

    for k1_ind = 1 : Nx
    for k2_ind = 1 : Ny

        eta_1 = kx_grid_2D(k1_ind, k2_ind);
        eta_2 = ky_grid_2D(k1_ind, k2_ind);

        q_1 = sqrt( (1 + dt * H44 * eta_2^2 + dt * G11 * eta_1^2) / (dt * H44) );
        q_2 = sqrt( (1 + dt * H44 * eta_1^2 + dt * G11 * eta_2^2) / (dt * H44) );
        q_3 = sqrt( (1 + dt * H44 * eta_1^2 + dt * H44 * eta_2^2) / (dt * G11) );

        q_mat_1(k1_ind,k2_ind) = q_1;
        q_mat_2(k1_ind,k2_ind) = q_2;
        q_mat_3(k1_ind,k2_ind) = q_3;

        scaling_1_pt = max(q_1 * z_axis(interface_index : film_index));
        scaling_2_pt = max(q_2 * z_axis(interface_index : film_index));
        scaling_3_pt = max(q_3 * z_axis(interface_index : film_index));

        mat = [ ...
            exp(q_1 * h_int - scaling_1_pt)   exp(-q_1 * h_int - scaling_1_pt); ...
            exp(q_1 * h_film - scaling_1_pt)  exp(-q_1 * h_film - scaling_1_pt) ];
        TimeCoeffs_inv_mat_1(k1_ind, k2_ind, :, :) = (1 / (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1))) * ...
                                                     [ mat(2,2)  -mat(1,2) ; ...
                                                       -mat(2,1) mat(1,1) ];

        mat = [ ...
            exp(q_2 * h_int - scaling_2_pt)   exp(-q_2 * h_int - scaling_2_pt); ...
            exp(q_2 * h_film - scaling_2_pt)  exp(-q_2 * h_film - scaling_2_pt) ];
        TimeCoeffs_inv_mat_2(k1_ind, k2_ind, :, :) = (1 / (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1))) * ...
                                                     [ mat(2,2)  -mat(1,2) ; ...
                                                       -mat(2,1) mat(1,1) ];

        mat = [ ...
            exp(q_3 * h_int - scaling_3_pt)   exp(-q_3 * h_int - scaling_3_pt); ...
            exp(q_3 * h_film - scaling_3_pt)  exp(-q_3 * h_film - scaling_3_pt) ];
        TimeCoeffs_inv_mat_3(k1_ind, k2_ind, :, :) = (1 / (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1))) * ...
                                                     [ mat(2,2)  -mat(1,2) ; ...
                                                       -mat(2,1) mat(1,1) ];

        scaling_1(k1_ind,k2_ind) = scaling_1_pt;
        scaling_2(k1_ind,k2_ind) = scaling_2_pt;
        scaling_3(k1_ind,k2_ind) = scaling_3_pt;

    end
    end
end