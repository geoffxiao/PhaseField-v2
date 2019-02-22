function [Thermal_Noise_1, Thermal_Noise_2, Thermal_Noise_3] = ThermalNoise(Constants)
    
    ThermConst = Constants.ThermConst;
    Nx = Constants.Nx;
    Ny = Constants.Ny;
    Nz = Constants.Nz;

    if(ThermConst ~= 0)
    % Thermal Noise
        Thermal_Noise_Sigma = sqrt(2 * ThermConst);

        Thermal_Noise_1 = normrnd(0,Thermal_Noise_Sigma,Nx,Ny,Nz);
        Thermal_Noise_2 = normrnd(0,Thermal_Noise_Sigma,Nx,Ny,Nz);
        Thermal_Noise_3 = normrnd(0,Thermal_Noise_Sigma,Nx,Ny,Nz);
    else
        Thermal_Noise_1 = 0;
        Thermal_Noise_2 = 0;
        Thermal_Noise_3 = 0;
	end
end