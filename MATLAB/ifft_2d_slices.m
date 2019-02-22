function out = ifft_2d_slices( mat )

    out = mat;
    [Nx, Ny, Nz] = size(mat);
    temp = zeros(Nx, Ny);
    for i = 1 : Nz
        temp(:) = mat(:,:,i);
        out(:,:,i) = real(ifft2(temp));
    end

end