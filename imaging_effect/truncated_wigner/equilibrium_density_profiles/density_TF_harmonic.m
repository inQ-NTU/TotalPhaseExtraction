function dens_profile = density_TF_harmonic(N_atoms, omegaZ_si, omegaT_si, grid_si)
% Calculate density profile of harmonic oscillator using the Thomas-Fermi 
% approximation, yielding an inverted parabola profile.
% Note, assumes oscillator to be centered on z=0

hbar_si     = 1.054571726e-34;  % reduced Planck's constant in SI-units
as_si       = 5.2e-9;           % scattering length in m
m_si        = 87*1.6605402e-27; % Rb87 mass

% calculate the coupling constant
g1d_si      = 2*hbar_si*as_si*omegaT_si;

% calculate central density
mu_TF_si    = 1/2*(m_si*(3/2*N_atoms*g1d_si*omegaZ_si)^2)^(1/3);
n0_si       = mu_TF_si/g1d_si;

% calculate Thomas-Fermi radius
R_TF        = sqrt(2*mu_TF_si/(m_si*omegaZ_si^2));

% density profile
dens_profile= n0_si*(1 - grid_si.^2/R_TF^2).*heaviside(R_TF - abs(grid));

end