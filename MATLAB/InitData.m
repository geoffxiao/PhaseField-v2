a = Ny / 4;

%% Init
P1 = zeros(Nx, Ny, Nz);
P2 = zeros(Nx, Ny, Nz);
P3 = zeros(Nx, Ny, Nz);

P1(:,a : 3 * a,:) = abs(rand(Nx,2 * a + 1,Nz));
P1(:,1 : a,:) = -abs(rand(Nx,a,Nz));
P1(:,3 * a : end,:) = -abs(rand(Nx,a + 1,Nz));

P2 = abs(rand(Nx,Ny,Nz));
P3 = abs(rand(Nx,Ny,Nz));

P1 = P1 * 1e-3;
P2 = P2 * 1e-3;
P3 = P3 * 1e-3;
