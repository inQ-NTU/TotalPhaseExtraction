function [cosine_Coeffs, sineCoeffs] = fourier_decomposition(signal, z_grid, max_modes)
    %extract condensate length and discretization length from the input grid
    L = abs(z_grid(end) - z_grid(1));
    dz = abs(z_grid(2) - z_grid(1));
    
    %fundamental wavevector
    k1 = 2*pi/L;

    %Setting the zero mode to zero
    %phase_profile = phase_profile - sum(phase_profile).*dz;

    %Finding cosineCoeffs
    for n = 1:max_modes
        coeff1 = 0;
        coeff2 = 0;
        for j = 1:length(z_grid)
            coeff1 = coeff1+signal(j)*cos(k1*n*z_grid(j));
            coeff2 = coeff2+signal(j)*sin(k1*n*z_grid(j));
        end
        cosine_Coeffs(n) = coeff1*(2*dz/L);
        sineCoeffs(n) = coeff2*(2*dz/L);
    end
end