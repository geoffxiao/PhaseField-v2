function [P1, P2, P3] = InitP(in_film, LOAD, Nx, Ny, Nz)

    % Initial Conditions
    if( LOAD )
        load('init.mat','P1','P2','P3');
        [n1, n2, n3] = size(P1);
                
        [yq,xq,zq] = meshgrid(linspace(1,n1,Ny),linspace(1,n2,Nx),linspace(1,n3,Nz));
        [y,x,z] = meshgrid(1:n1,1:n2,1:n3);
        
        P1 = interp3(y,x,z,P1,yq,xq,zq);
        P2 = interp3(y,x,z,P2,yq,xq,zq);
        P3 = interp3(y,x,z,P3,yq,xq,zq);
      
    else
        [P1, P2, P3] = RandomP(Nx, Ny, Nz);
    end

    % % Thin film simulation
    P1 = P1 .* in_film;
    P2 = P2 .* in_film;
    P3 = P3 .* in_film;


% a = Ny / 4;
% 
% % Init
% P1 = zeros(Nx, Ny, Nz);
% P2 = zeros(Nx, Ny, Nz);
% P3 = zeros(Nx, Ny, Nz);
% 
% P1(:,a : 3 * a,:) = abs(rand(Nx,2 * a + 1,Nz));
% P1(:,1 : a,:) = -abs(rand(Nx,a,Nz));
% P1(:,3 * a : end,:) = -abs(rand(Nx,a + 1,Nz));
% 
% P2 = abs(rand(Nx,Ny,Nz));
% P3 = abs(rand(Nx,Ny,Nz));
% 
% P1 = P1 * 1e-3;
% P2 = P2 * 1e-3;
% P3 = P3 * 1e-3;
% 
% P1 = P1 .* in_film;
% P2 = P2 .* in_film;
% P3 = P3 .* in_film;
    
end