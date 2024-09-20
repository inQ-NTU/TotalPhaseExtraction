function modes_Coeffs = cosine_decomposition(phase_profile, z_grid, max_modes)
    %extract condensate length and discretization length from the input grid
    L = abs(z_grid(end) - z_grid(1));
    dz = abs(z_grid(2) - z_grid(1));
    
    %fundamental wavevector
    k1 = pi/L;

    %Setting the zero mode to zero
    %phase_profile = phase_profile - sum(phase_profile).*dz;

    %Finding cosineCoeffs
    for n = 1:max_modes
        coeff = 0;
        for j = 1:length(z_grid)
            coeff = coeff+phase_profile(j)*cos(k1*n*z_grid(j));
        end
        modes_Coeffs(n) = coeff*(2*dz/L);
    end
end