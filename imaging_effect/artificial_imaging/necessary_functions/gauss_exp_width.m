function sigma_n = gauss_exp_width(omega_perp, n, t_tof)
% width of the expanding transverse wfnk assuming a Gaussian shape and that
% all initial kinetic and interaction energy goes into the expansion

hbar = 1.05457148e-34;          % reduced plancks constant
as = 98.98*52.917720859e-12;    % s-wave scattering length
m = 1.443160e-25;   % Rb87 mass

sigma_0 = sqrt(hbar/(m*omega_perp)); %harmonic oscillator length

sigma_n = sigma_0 * (1 + 2*as*n).^(1/4) * sqrt(1 + (omega_perp*t_tof).^2);

end