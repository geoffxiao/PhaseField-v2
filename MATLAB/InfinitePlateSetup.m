function [strain_bc_mats_inv, strain_bc_mat_inv_korigin, ...
    strain_sol_2Dk_1, strain_sol_2Dk_2, strain_sol_2Dk_3, ...
    strain_sol_2Dk_4, strain_sol_2Dk_5, strain_sol_2Dk_6, ...
    strain_sol_2Dk_d3_1, strain_sol_2Dk_d3_2, strain_sol_2Dk_d3_3, ...
    strain_sol_2Dk_d3_4, strain_sol_2Dk_d3_5, strain_sol_2Dk_d3_6, ...
    strain_bc_mats_L_inv_lu, strain_bc_mats_U_lu, strain_bc_mats_P_lu] = InfinitePlateSetup(Constants, kx_grid_2D, ky_grid_2D)
    
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    Nz = Constants.Nz;
    C = Constants.C;
    kx = Constants.kx;
    ky = Constants.ky;
    h_film = Constants.h_film;
    h_sub = Constants.h_sub;
    z_axis = Constants.z_axis;
                       
    U = zeros(3,3);
    m = [0, 0, 1];
    for i_ind = 1 : 3
    for k_ind = 1 : 3
        for j_ind = 1 : 3
        for l_ind = 1 : 3
            U(i_ind,k_ind) = U(i_ind,k_ind) + C(i_ind,j_ind,k_ind,l_ind) * m(j_ind) * m(l_ind);
        end
        end
    end
    end

    R = zeros(Nx,Ny,3,3); W = zeros(numel(kx),numel(ky),3,3);
    eigenvec_mat = zeros(Nx,Ny,6,6); % the eigenvalues and eigenvectors found at each k vector pt
    eigenval_mat = zeros(Nx,Ny,6); % the eigenvalues and eigenvectors found at each k vector pt
    mat_mat = zeros(numel(kx),numel(ky),6,6);

    for k1_ind = 1 : numel(kx)
    for k2_ind = 1 : numel(ky)

        if( k1_ind ~= 1 || k2_ind ~= 1 ) % Origin, different solution, be mindful of dividing by k vector as well!

            n = [ kx_grid_2D(k1_ind,k2_ind), ky_grid_2D(k1_ind,k2_ind), 0 ] ./ sqrt( kx_grid_2D(k1_ind,k2_ind)^2 + ky_grid_2D(k1_ind,k2_ind)^2 ) ;

            for i_ind = 1 : 3
            for k_ind = 1 : 3
                for j_ind = 1 : 3
                for l_ind = 1 : 3
                    R(k1_ind,k2_ind,i_ind,k_ind) = R(k1_ind,k2_ind,i_ind,k_ind) + C(i_ind,j_ind,k_ind,l_ind) * n(j_ind) * m(l_ind);
                    W(k1_ind,k2_ind,i_ind,k_ind) = W(k1_ind,k2_ind,i_ind,k_ind) + C(i_ind,j_ind,k_ind,l_ind) * n(j_ind) * n(l_ind);
                end
                end
            end
            end

            R_mat = squeeze(R(k1_ind,k2_ind,:,:));
            W_mat = squeeze(W(k1_ind,k2_ind,:,:));


            N1 = -inv(U)*R_mat'; N2 = inv(U); N3 = R_mat*inv(U)*R_mat'-W_mat;
            mat = [N1, N2; N3, N1'];

            [eigenvec, eigenval] = eig(mat,'vector');

            mat_mat(k1_ind,k2_ind,:,:) = mat;
            eigenvec_mat(k1_ind,k2_ind,:,:) = eigenvec;
            eigenval_mat(k1_ind,k2_ind,:) = eigenval;

        end
    end
    end

    %% BC matrix solutions
    bc_mats = zeros(Nx,Ny,6,6);
    strain_bc_mats_inv = zeros(Nx,Ny,6,6);

    % LU Decomp
    strain_bc_mats_U_lu = zeros(Nx,Ny,6,6);
    strain_bc_mats_L_lu = zeros(Nx,Ny,6,6);
    strain_bc_mats_P_lu = zeros(Nx,Ny,6,6);
    strain_bc_mats_L_inv_lu = zeros(Nx,Ny,6,6);
    
    for k1_ind = 1 : Nx
    for k2_ind = 1 : Ny
        if( k1_ind ~= 1 || k2_ind ~= 1 ) % Because we divide by K^2 to make matrices better
            K = sqrt( kx_grid_2D(k1_ind,k2_ind)^2 + ky_grid_2D(k1_ind,k2_ind)^2 );
            bc_mat = zeros(6,6);
            % [ d1 * eigenvec(1,1) + d2 * eigenvec(1,2) + 
            %   d1 * eigenvec(2,1) + d2 * eigenvec(2,2) + 
            %   d1 * eigenvec(3,1) + d2 * eigenvec(3,2) + 
            %   d1 * eigenvec(4,1) + d2 * eigenvec(4,2) + 
            for m = 1 : 6 % the eigenvalue -> column of eigenvec_mat ~ eigenvector
                for k_ind = 1 : 3 % eigenvector a, scale by K 
                    bc_mat(k_ind,m) = eigenvec_mat(k1_ind,k2_ind,k_ind,m) * exp(eigenval_mat(k1_ind,k2_ind,m)*h_sub*1i*K);
                end
                % rescale by K^-1 to avoid conditioning errors... need to rescale
                % the given value as well!
                for k_ind = 4 : 6 % eigenvector b
                    bc_mat(k_ind,m) = eigenvec_mat(k1_ind,k2_ind,k_ind,m) * exp(eigenval_mat(k1_ind,k2_ind,m)*h_film*1i*K) * 1i * K / K^2;        
                end
            end
            bc_mats(k1_ind,k2_ind,:,:) = bc_mat;
            % strain_bc_mats_inv(k1_ind,k2_ind,:,:) = inv(bc_mat);

            [L_lu,U_lu,P_lu] = lu(bc_mat);
            strain_bc_mats_L_lu(k1_ind,k2_ind,:,:) = L_lu;
            strain_bc_mats_L_inv_lu(k1_ind,k2_ind,:,:) = inv(L_lu);
            strain_bc_mats_U_lu(k1_ind,k2_ind,:,:) = U_lu;
            strain_bc_mats_P_lu(k1_ind,k2_ind,:,:) = P_lu;
                
        end
    end
    end

    strain_bc_mat_inv_korigin = inv([ h_sub        1   0           0   0           0 ; ...
                                      0            0   h_sub       1   0           0 ; ...
                                      0            0   0           0   h_sub       1 ; ...
                                      C(1,3,1,3)   0   C(1,3,2,3)  0   C(1,3,3,3)  0 ; ...
                                      C(2,3,1,3)   0   C(2,3,2,3)  0   C(2,3,3,3)  0 ; ...
                                      C(3,3,1,3)   0   C(3,3,2,3)  0   C(3,3,3,3)  0 ]);
 
    strain_sol_2Dk_1 = zeros(Nx, Ny, Nz, 3);                              
    strain_sol_2Dk_2 = zeros(Nx, Ny, Nz, 3);                              
    strain_sol_2Dk_3 = zeros(Nx, Ny, Nz, 3);                              
    strain_sol_2Dk_4 = zeros(Nx, Ny, Nz, 3);                              
    strain_sol_2Dk_5 = zeros(Nx, Ny, Nz, 3);                              
    strain_sol_2Dk_6 = zeros(Nx, Ny, Nz, 3);                         

    strain_sol_2Dk_d3_1 = zeros(Nx, Ny, Nz, 3);    
    strain_sol_2Dk_d3_2 = zeros(Nx, Ny, Nz, 3);    
    strain_sol_2Dk_d3_3 = zeros(Nx, Ny, Nz, 3);    
    strain_sol_2Dk_d3_4 = zeros(Nx, Ny, Nz, 3);    
    strain_sol_2Dk_d3_5 = zeros(Nx, Ny, Nz, 3);    
    strain_sol_2Dk_d3_6 = zeros(Nx, Ny, Nz, 3);    

    p_1 = squeeze(eigenval_mat(:,:,1));
    p_2 = squeeze(eigenval_mat(:,:,2));
    p_3 = squeeze(eigenval_mat(:,:,3));
    p_4 = squeeze(eigenval_mat(:,:,4));
    p_5 = squeeze(eigenval_mat(:,:,5));
    p_6 = squeeze(eigenval_mat(:,:,6));

    for l = 1 : 3
        for z_loop = 1 : numel(z_axis);
            K = sqrt( kx_grid_2D.^2 + ky_grid_2D.^2 );

            a_vec_1 = squeeze(eigenvec_mat(:,:,l,1));
            a_vec_2 = squeeze(eigenvec_mat(:,:,l,2));
            a_vec_3 = squeeze(eigenvec_mat(:,:,l,3));
            a_vec_4 = squeeze(eigenvec_mat(:,:,l,4));
            a_vec_5 = squeeze(eigenvec_mat(:,:,l,5));
            a_vec_6 = squeeze(eigenvec_mat(:,:,l,6));

            strain_sol_2Dk_1(:,:,z_loop,l) = a_vec_1 .* exp(p_1.*z_axis(z_loop).*1i.*K);
            strain_sol_2Dk_2(:,:,z_loop,l) = a_vec_2 .* exp(p_2.*z_axis(z_loop).*1i.*K);
            strain_sol_2Dk_3(:,:,z_loop,l) = a_vec_3 .* exp(p_3.*z_axis(z_loop).*1i.*K);
            strain_sol_2Dk_4(:,:,z_loop,l) = a_vec_4 .* exp(p_4.*z_axis(z_loop).*1i.*K);
            strain_sol_2Dk_5(:,:,z_loop,l) = a_vec_5 .* exp(p_5.*z_axis(z_loop).*1i.*K);
            strain_sol_2Dk_6(:,:,z_loop,l) = a_vec_6 .* exp(p_6.*z_axis(z_loop).*1i.*K);

            strain_sol_2Dk_d3_1(:,:,z_loop,l) = a_vec_1 .* exp(p_1.*z_axis(z_loop).*1i.*K) .* p_1 .* 1i .* K;
            strain_sol_2Dk_d3_2(:,:,z_loop,l) = a_vec_2 .* exp(p_2.*z_axis(z_loop).*1i.*K) .* p_2 .* 1i .* K;
            strain_sol_2Dk_d3_3(:,:,z_loop,l) = a_vec_3 .* exp(p_3.*z_axis(z_loop).*1i.*K) .* p_3 .* 1i .* K;
            strain_sol_2Dk_d3_4(:,:,z_loop,l) = a_vec_4 .* exp(p_4.*z_axis(z_loop).*1i.*K) .* p_4 .* 1i .* K;
            strain_sol_2Dk_d3_5(:,:,z_loop,l) = a_vec_5 .* exp(p_5.*z_axis(z_loop).*1i.*K) .* p_5 .* 1i .* K;
            strain_sol_2Dk_d3_6(:,:,z_loop,l) = a_vec_6 .* exp(p_6.*z_axis(z_loop).*1i.*K) .* p_6 .* 1i .* K; 
        end
    end
end