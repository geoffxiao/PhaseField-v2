function [f1_elastic, f2_elastic, f3_elastic, u_1, u_2, u_3, ...
    TotalStrain_11, TotalStrain_12, TotalStrain_13, ...
    TotalStrain_22, TotalStrain_23, TotalStrain_33] = ...
            CalcElasticEnergy(P1, P2, P3, Constants, ...
            Green_11, Green_12, Green_13, ...
            Green_21, Green_22, Green_23, ...
            Green_31, Green_32, Green_33, ...
            strain_bc_mats_inv, strain_bc_mat_inv_korigin, ...
            kx_grid_3D, ky_grid_3D, kz_grid_3D, z_axis, ...
            strain_sol_2Dk_1, strain_sol_2Dk_2, strain_sol_2Dk_3, ...
            strain_sol_2Dk_4, strain_sol_2Dk_5, strain_sol_2Dk_6, ...
            strain_sol_2Dk_d3_1, strain_sol_2Dk_d3_2, strain_sol_2Dk_d3_3, ...
            strain_sol_2Dk_d3_4, strain_sol_2Dk_d3_5, strain_sol_2Dk_d3_6, ...
            strain_bc_mats_L_inv_lu, strain_bc_mats_U_lu, strain_bc_mats_P_lu)

    C11 = Constants.C11;
    C12 = Constants.C12;
    C44 = Constants.C44;
    Q11 = Constants.Q11;
    Q12 = Constants.Q12;
    Q44 = Constants.Q44;   
    q11 = Constants.q11;
    q12 = Constants.q12;
    q44 = Constants.q44;
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    Nz = Constants.Nz;
    dx = Constants.dx;
    dy = Constants.dy;
    dz = Constants.dz;
    film_index = Constants.film_index;
    ELASTIC = Constants.ELASTIC;
    HET_ELASTIC_RELAX = Constants.HET_ELASTIC_RELAX;
    
    kx_grid_2D = kx_grid_3D(:,:,1);
    ky_grid_2D = ky_grid_3D(:,:,1);    
    
    if( ELASTIC )

        if(HET_ELASTIC_RELAX) % Heterogenous relaxation

            % Eigenstrains
            Eigenstrain_11 = Q11 * P1.^2 + Q12 .* (P2.^2 + P3.^2);
            Eigenstrain_22 = Q11 * P2.^2 + Q12 .* (P1.^2 + P3.^2);
            Eigenstrain_33 = Q11 * P3.^2 + Q12 .* (P1.^2 + P2.^2);
            Eigenstrain_23 = Q44 * P2 .* P3;
            Eigenstrain_13 = Q44 * P1 .* P3;
            Eigenstrain_12 = Q44 * P1 .* P2; 

            Eigenstrain_11_k = fftn(Eigenstrain_11);
            Eigenstrain_22_k = fftn(Eigenstrain_22);
            Eigenstrain_33_k = fftn(Eigenstrain_33);
            Eigenstrain_23_k = fftn(Eigenstrain_23);
            Eigenstrain_13_k = fftn(Eigenstrain_13);
            Eigenstrain_12_k = fftn(Eigenstrain_12);
            
            % Calculate the displacement field via Khachutaryan method
            % For each k vector calculate displacement field
            u_1_A_k = - C11.*Eigenstrain_11_k.*kx_grid_3D.*Green_11.*1i - ...
                        C12.*Eigenstrain_22_k.*kx_grid_3D.*Green_11.*1i - ...
                        C12.*Eigenstrain_33_k.*kx_grid_3D.*Green_11.*1i - ...
                        C44.*Eigenstrain_12_k.*kx_grid_3D.*Green_12.*2i - ...
                        C44.*Eigenstrain_13_k.*kx_grid_3D.*Green_13.*2i - ...
                        C12.*Eigenstrain_11_k.*ky_grid_3D.*Green_12.*1i - ...
                        C11.*Eigenstrain_22_k.*ky_grid_3D.*Green_12.*1i - ...
                        C12.*Eigenstrain_33_k.*ky_grid_3D.*Green_12.*1i - ...
                        C44.*Eigenstrain_12_k.*ky_grid_3D.*Green_11.*2i - ...
                        C44.*Eigenstrain_23_k.*ky_grid_3D.*Green_13.*2i - ...
                        C12.*Eigenstrain_11_k.*kz_grid_3D.*Green_13.*1i - ...
                        C12.*Eigenstrain_22_k.*kz_grid_3D.*Green_13.*1i - ...
                        C11.*Eigenstrain_33_k.*kz_grid_3D.*Green_13.*1i - ...
                        C44.*Eigenstrain_13_k.*kz_grid_3D.*Green_11.*2i - ...
                        C44.*Eigenstrain_23_k.*kz_grid_3D.*Green_12.*2i;
            u_2_A_k = - C11.*Eigenstrain_11_k.*kx_grid_3D.*Green_21.*1i - ...
                        C12.*Eigenstrain_22_k.*kx_grid_3D.*Green_21.*1i - ...
                        C12.*Eigenstrain_33_k.*kx_grid_3D.*Green_21.*1i - ...
                        C44.*Eigenstrain_12_k.*kx_grid_3D.*Green_22.*2i - ...
                        C44.*Eigenstrain_13_k.*kx_grid_3D.*Green_23.*2i - ...
                        C12.*Eigenstrain_11_k.*ky_grid_3D.*Green_22.*1i - ...
                        C11.*Eigenstrain_22_k.*ky_grid_3D.*Green_22.*1i - ...
                        C12.*Eigenstrain_33_k.*ky_grid_3D.*Green_22.*1i - ...
                        C44.*Eigenstrain_12_k.*ky_grid_3D.*Green_21.*2i - ...
                        C44.*Eigenstrain_23_k.*ky_grid_3D.*Green_23.*2i - ...
                        C12.*Eigenstrain_11_k.*kz_grid_3D.*Green_23.*1i - ...
                        C12.*Eigenstrain_22_k.*kz_grid_3D.*Green_23.*1i - ...
                        C11.*Eigenstrain_33_k.*kz_grid_3D.*Green_23.*1i - ...
                        C44.*Eigenstrain_13_k.*kz_grid_3D.*Green_21.*2i - ...
                        C44.*Eigenstrain_23_k.*kz_grid_3D.*Green_22.*2i;
            u_3_A_k = - C11.*Eigenstrain_11_k.*kx_grid_3D.*Green_31.*1i - ...
                        C12.*Eigenstrain_22_k.*kx_grid_3D.*Green_31.*1i - ...
                        C12.*Eigenstrain_33_k.*kx_grid_3D.*Green_31.*1i - ...
                        C44.*Eigenstrain_12_k.*kx_grid_3D.*Green_32.*2i - ...
                        C44.*Eigenstrain_13_k.*kx_grid_3D.*Green_33.*2i - ...
                        C12.*Eigenstrain_11_k.*ky_grid_3D.*Green_32.*1i - ...
                        C11.*Eigenstrain_22_k.*ky_grid_3D.*Green_32.*1i - ...
                        C12.*Eigenstrain_33_k.*ky_grid_3D.*Green_32.*1i - ...
                        C44.*Eigenstrain_12_k.*ky_grid_3D.*Green_31.*2i - ...
                        C44.*Eigenstrain_23_k.*ky_grid_3D.*Green_33.*2i - ...
                        C12.*Eigenstrain_11_k.*kz_grid_3D.*Green_33.*1i - ...
                        C12.*Eigenstrain_22_k.*kz_grid_3D.*Green_33.*1i - ...
                        C11.*Eigenstrain_33_k.*kz_grid_3D.*Green_33.*1i - ...
                        C44.*Eigenstrain_13_k.*kz_grid_3D.*Green_31.*2i - ...
                        C44.*Eigenstrain_23_k.*kz_grid_3D.*Green_32.*2i;
        
            zero_index = kx_grid_3D == 0 & ky_grid_3D == 0 & ky_grid_3D == 0;
 
            % ifftn u_i
            % DC component = zero to allow ifft
            u_1_A_k(zero_index) = 0;
            u_2_A_k(zero_index) = 0;
            u_3_A_k(zero_index) = 0;
            
            e_11_A_k = 1i .* kx_grid_3D .* u_1_A_k;
            e_22_A_k = 1i .* ky_grid_3D .* u_2_A_k;
            e_33_A_k = 1i .* kz_grid_3D .* u_3_A_k;
            e_12_A_k = 0.5 * ( 1i .* kx_grid_3D .* u_2_A_k + 1i .* ky_grid_3D .* u_1_A_k );
            e_13_A_k = 0.5 * ( 1i .* kx_grid_3D .* u_3_A_k + 1i .* kz_grid_3D .* u_1_A_k );
            e_23_A_k = 0.5 * ( 1i .* kz_grid_3D .* u_2_A_k + 1i .* ky_grid_3D .* u_3_A_k );

            % real space ifftn
            u_1_A = real(ifftn(u_1_A_k));
            u_2_A = real(ifftn(u_2_A_k));
            u_3_A = real(ifftn(u_3_A_k));
            
            e_11_A = real(ifftn(e_11_A_k));
            e_22_A = real(ifftn(e_22_A_k));
            e_33_A = real(ifftn(e_33_A_k));
            e_13_A = real(ifftn(e_13_A_k));
            e_23_A = real(ifftn(e_23_A_k));
            e_12_A = real(ifftn(e_12_A_k));

            %% Application of Boundary Conditions
            % 3D fft diff..
            bc_film_1_3Dk = C44.*(Eigenstrain_13_k - kx_grid_3D.*u_3_A_k.*1i) + C44.*(Eigenstrain_13_k - kz_grid_3D.*u_1_A_k.*1i);
            bc_film_2_3Dk = C44.*(Eigenstrain_23_k - ky_grid_3D.*u_3_A_k.*1i) + C44.*(Eigenstrain_23_k - kz_grid_3D.*u_2_A_k.*1i);
            bc_film_3_3Dk = C12.*(Eigenstrain_11_k - kx_grid_3D.*u_1_A_k.*1i) + C12.*(Eigenstrain_22_k - ky_grid_3D.*u_2_A_k.*1i) + C11.*(Eigenstrain_33_k - kz_grid_3D.*u_3_A_k.*1i);

            % real space bc
            % bc applied on the top surface (surface of the film) in real space
            bc_film_1 = real(ifftn(squeeze(bc_film_1_3Dk)));
            bc_film_2 = real(ifftn(squeeze(bc_film_2_3Dk)));
            bc_film_3 = real(ifftn(squeeze(bc_film_3_3Dk)));

            bc_film_1 = squeeze(bc_film_1(:,:,film_index));
            bc_film_2 = squeeze(bc_film_2(:,:,film_index));
            bc_film_3 = squeeze(bc_film_3(:,:,film_index));

            % Strain free substrate
            bc_sub_1 = -squeeze(u_1_A(:,:,1));
            bc_sub_2 = -squeeze(u_2_A(:,:,1));
            bc_sub_3 = -squeeze(u_3_A(:,:,1));

            %% fft2 the bc 
            bc_film_1_2Dk = fft2(bc_film_1);
            bc_film_2_2Dk = fft2(bc_film_2);
            bc_film_3_2Dk = fft2(bc_film_3);

            % Orders of magnitude larger... so let's make the matrix better to work with
            rescaling_mat = sqrt( kx_grid_2D.^2 + ky_grid_2D.^2 ).^2;
            rescaling_mat(1,1) = 1;
            bc_film_1_2Dk = bc_film_1_2Dk ./ rescaling_mat;
            bc_film_2_2Dk = bc_film_2_2Dk ./ rescaling_mat;
            bc_film_3_2Dk = bc_film_3_2Dk ./ rescaling_mat;

            bc_sub_1_2Dk = fft2(bc_sub_1);
            bc_sub_2_2Dk = fft2(bc_sub_2);
            bc_sub_3_2Dk = fft2(bc_sub_3);

            bc_given_k_space = cat(3,bc_sub_1_2Dk,bc_sub_2_2Dk,bc_sub_3_2Dk,bc_film_1_2Dk,bc_film_2_2Dk,bc_film_3_2Dk);
            d_sols = zeros(Nx,Ny,6);

            %% Apply the boundary conditions to get the solution to the diff. eq.
            temp1 = zeros(6, 6);
            temp2 = zeros(6, 6);
            temp3 = zeros(6, 6);
            temp4 = zeros(6, 1);
            y = zeros(6, 1);
            for k1 = 1 : Nx
            for k2 = 1 : Ny
            if(k1 ~= 1 || k2 ~= 1)
                % Solve the matrix problem
            %    temp1(:) = strain_bc_mats_inv(k1,k2,:,:);
            %    temp2(:) = bc_given_k_space(k1,k2,:);
            %    d_sols(k1,k2,:) = squeeze(strain_bc_mats_inv(k1,k2,:,:)) * squeeze(bc_given_k_space(k1,k2,:));
            %    d_sols(k1,k2,:) = temp1 * temp2;
            
                % 1 = L, 2 = U, 3 = P, 4 = bc_given
                temp1(:) = strain_bc_mats_L_inv_lu(k1,k2,:,:);
                temp2(:) = strain_bc_mats_U_lu(k1,k2,:,:);
                temp3(:) = strain_bc_mats_P_lu(k1,k2,:,:);
                temp4(:) = bc_given_k_space(k1,k2,:);
                y = temp1 * temp3 * temp4;
                d_sols(k1,k2,:) = temp2 \ y;
            end
            end
            end

            %% Construct the solution in the 2D k-space, kx,ky,z
            u_B_2Dk_mat = zeros(Nx,Ny,Nz,3);
            u_B_2Dk_d3_mat = zeros(Nx,Ny,Nz,3);

            q_1 = squeeze(d_sols(:,:,1));
            q_2 = squeeze(d_sols(:,:,2));
            q_3 = squeeze(d_sols(:,:,3));
            q_4 = squeeze(d_sols(:,:,4));
            q_5 = squeeze(d_sols(:,:,5));
            q_6 = squeeze(d_sols(:,:,6));

            q_1 = repmat(q_1,1,1,Nz);
            q_2 = repmat(q_2,1,1,Nz);
            q_3 = repmat(q_3,1,1,Nz);
            q_4 = repmat(q_4,1,1,Nz);
            q_5 = repmat(q_5,1,1,Nz);
            q_6 = repmat(q_6,1,1,Nz);

            s1 = zeros(Nx,Ny,Nz);
            s2 = zeros(Nx,Ny,Nz);
            s3 = zeros(Nx,Ny,Nz);
            s4 = zeros(Nx,Ny,Nz);
            s5 = zeros(Nx,Ny,Nz);
            s6 = zeros(Nx,Ny,Nz);

            s1_d3 = zeros(Nx,Ny,Nz);
            s2_d3 = zeros(Nx,Ny,Nz);
            s3_d3 = zeros(Nx,Ny,Nz);
            s4_d3 = zeros(Nx,Ny,Nz);
            s5_d3 = zeros(Nx,Ny,Nz);
            s6_d3 = zeros(Nx,Ny,Nz);

            for l = 1 : 3
                    s1(:) = strain_sol_2Dk_1(:,:,:,l);
                    s2(:) = strain_sol_2Dk_2(:,:,:,l);
                    s3(:) = strain_sol_2Dk_3(:,:,:,l);
                    s4(:) = strain_sol_2Dk_4(:,:,:,l);
                    s5(:) = strain_sol_2Dk_5(:,:,:,l);
                    s6(:) = strain_sol_2Dk_6(:,:,:,l);

                    s1_d3(:) = strain_sol_2Dk_d3_1(:,:,:,l);
                    s2_d3(:) = strain_sol_2Dk_d3_2(:,:,:,l);
                    s3_d3(:) = strain_sol_2Dk_d3_3(:,:,:,l);
                    s4_d3(:) = strain_sol_2Dk_d3_4(:,:,:,l);
                    s5_d3(:) = strain_sol_2Dk_d3_5(:,:,:,l);
                    s6_d3(:) = strain_sol_2Dk_d3_6(:,:,:,l);

                    u_B_2Dk_mat(:,:,:,l) = q_1 .* s1 + ...
                                                q_2 .* s2 + ...
                                                q_3 .* s3 + ...
                                                q_4 .* s4 + ...
                                                q_5 .* s5 + ...
                                                q_6 .* s6;

                    u_B_2Dk_d3_mat(:,:,:,l) = q_1 .* s1_d3 + ...
                                                   q_2 .* s2_d3 + ...
                                                   q_3 .* s3_d3 + ...
                                                   q_4 .* s4_d3 + ...
                                                   q_5 .* s5_d3 + ...
                                                   q_6 .* s6_d3;
            end

            % 2D origin point...
            d_sols(1,1,:) = strain_bc_mat_inv_korigin * squeeze(bc_given_k_space(1,1,:));

            q = squeeze(d_sols(1,1,:));
            u_B_2Dk_mat(1,1,:,1) = q(1) * z_axis + q(2);
            u_B_2Dk_mat(1,1,:,2) = q(3) * z_axis + q(4);
            u_B_2Dk_mat(1,1,:,3) = q(5) * z_axis + q(6);
            u_B_2Dk_d3_mat(1,1,:,1) = q(1);
            u_B_2Dk_d3_mat(1,1,:,2) = q(3);
            u_B_2Dk_d3_mat(1,1,:,3) = q(5);

            u_1_B_2Dk = squeeze(u_B_2Dk_mat(:,:,:,1));
            u_2_B_2Dk = squeeze(u_B_2Dk_mat(:,:,:,2));
            u_3_B_2Dk = squeeze(u_B_2Dk_mat(:,:,:,3));
            u_1_B_2Dk_d3 = squeeze(u_B_2Dk_d3_mat(:,:,:,1));
            u_2_B_2Dk_d3 = squeeze(u_B_2Dk_d3_mat(:,:,:,2));
            u_3_B_2Dk_d3 = squeeze(u_B_2Dk_d3_mat(:,:,:,3));

            %% Calculate the strains, e_ij = 0.5 * (u_i,j + u_j,i)
            e_11_B_2Dk = u_1_B_2Dk .* 1i .* kx_grid_3D;
            e_22_B_2Dk = u_2_B_2Dk .* 1i .* ky_grid_3D;
            e_33_B_2Dk = u_3_B_2Dk_d3;
            e_23_B_2Dk = 0.5 * (u_2_B_2Dk_d3 + 1i .* ky_grid_3D .* u_3_B_2Dk);
            e_13_B_2Dk = 0.5 * (u_1_B_2Dk_d3 + 1i .* kx_grid_3D .* u_3_B_2Dk);
            e_12_B_2Dk = 0.5 * (1i .* ky_grid_3D .* u_1_B_2Dk + 1i .* kx_grid_3D .* u_2_B_2Dk);

            %% Final outputs
            u_1_B = ifft_2d_slices(squeeze(u_B_2Dk_mat(:,:,:,1)));
            u_2_B = ifft_2d_slices(squeeze(u_B_2Dk_mat(:,:,:,2)));
            u_3_B = ifft_2d_slices(squeeze(u_B_2Dk_mat(:,:,:,3)));

            e_11_B = ifft_2d_slices(e_11_B_2Dk);
            e_22_B = ifft_2d_slices(e_22_B_2Dk);
            e_33_B = ifft_2d_slices(e_33_B_2Dk);
            e_23_B = ifft_2d_slices(e_23_B_2Dk);
            e_13_B = ifft_2d_slices(e_13_B_2Dk);
            e_12_B = ifft_2d_slices(e_12_B_2Dk);


        else % Don't do heterogenous strain relaxation
            u_1_A = 0;
            u_2_A = 0;
            u_3_A = 0;

            u_1_B = 0;
            u_2_B = 0;
            u_3_B = 0;    

            e_11_A = 0;
            e_22_A = 0;
            e_33_A = 0;
            e_12_A = 0;    
            e_13_A = 0;
            e_23_A = 0;

            e_11_B = 0;
            e_22_B = 0;
            e_33_B = 0;
            e_12_B = 0;    
            e_13_B = 0;
            e_23_B = 0;
        end

        u_1 = u_1_A + u_1_B; u_2 = u_2_A + u_2_B; u_3 = u_3_A + u_3_B;
        
        e_11_het = e_11_A + e_11_B;
        e_22_het = e_22_A + e_22_B;
        e_33_het = e_33_A + e_33_B;
        e_12_het = e_12_A + e_12_B;
        e_13_het = e_13_A + e_13_B;
        e_23_het = e_23_A + e_23_B;
        
        % rescale so volume integrao of het strain = 0
        e_11_het = e_11_het - ( trapz(trapz(trapz(e_11_het))) / ( Nx * Ny * Nz ));
        e_22_het = e_22_het - ( trapz(trapz(trapz(e_22_het))) / ( Nx * Ny * Nz ));
        e_33_het = e_33_het - ( trapz(trapz(trapz(e_33_het))) / ( Nx * Ny * Nz ));
        e_12_het = e_12_het - ( trapz(trapz(trapz(e_12_het))) / ( Nx * Ny * Nz ));
        e_13_het = e_13_het - ( trapz(trapz(trapz(e_13_het))) / ( Nx * Ny * Nz ));
        e_23_het = e_23_het - ( trapz(trapz(trapz(e_23_het))) / ( Nx * Ny * Nz ));
        
        TotalStrain_11 = e_11_het + TotalStrain_homo_11(P1,P2,P3,Constants); 
        TotalStrain_22 = e_22_het + TotalStrain_homo_22(P1,P2,P3,Constants); 
        TotalStrain_33 = e_33_het + TotalStrain_homo_33(P1,P2,P3,Constants);
        TotalStrain_12 = e_12_het + TotalStrain_homo_12(P1,P2,P3,Constants);
        TotalStrain_23 = e_23_het + TotalStrain_homo_23(P1,P2,P3,Constants); 
        TotalStrain_13 = e_13_het + TotalStrain_homo_13(P1,P2,P3,Constants); 
        
        
        %% Elastic Energy
        b11 = 0.5*C11*(Q11^2 + 2*Q12^2) + C12*Q12*(2*Q11 + Q12);
        b12 = C11*Q12*(2*Q11 + Q12) + C12*(Q11^2 + 3*Q12^2 + 2*Q11*Q12) + 2*C44*Q44^2;

        f1_elastic = -( (q11.*TotalStrain_11 + q12.*TotalStrain_22 + q12.*TotalStrain_33) .* (2.*P1) + 2.*q44.*(TotalStrain_12 .* P2 + TotalStrain_13 .* P3) ) + ( 4 .* b11 .* P1.^2 + 2 .* b12 .* ( P2.^2 + P3.^2 ) ) .* P1;
        f2_elastic = -( (q11.*TotalStrain_22 + q12.*TotalStrain_11 + q12.*TotalStrain_33) .* (2.*P2) + 2.*q44.*(TotalStrain_12 .* P1 + TotalStrain_13 .* P3) ) + ( 4 .* b11 .* P2.^2 + 2 .* b12 .* ( P3.^2 + P1.^2 ) ) .* P2;
        f3_elastic = -( (q11.*TotalStrain_33 + q12.*TotalStrain_11 + q12.*TotalStrain_22) .* (2.*P3) + 2.*q44.*(TotalStrain_23 .* P2 + TotalStrain_13 .* P1) ) + ( 4 .* b11 .* P3.^2 + 2 .* b12 .* ( P1.^2 + P2.^2 ) ) .* P3;

    else
        
        f1_elastic = 0; f2_elastic = 0; f3_elastic = 0;
        u_1 = 0; u_2 = 0; u_3 = 0;
        
    end
    
end