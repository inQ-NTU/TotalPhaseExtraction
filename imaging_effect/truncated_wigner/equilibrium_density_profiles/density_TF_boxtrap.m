function dens_profile = density_TF_boxtrap(N_atoms, boxlength, psf_DMD, grid_si)
% Calculate density profile of boxtrap using the Thomas-Fermi 
% approximation, yielding an inverted parabola profile.
% Note, assumes box to be centered on z=0

hbar_si     = 1.054571726e-34;  % reduced Planck's constant in SI-units
as_si       = 5.2e-9;           % scattering length in m
m_si        = 87*1.6605402e-27; % Rb87 mass

% calculate central density
n0_si       = N_atoms/boxlength;

% density profile
dens_profile = n0_si*(heaviside(grid_si + boxlength/2)- heaviside(grid_si - boxlength/2));

% account for PSF of DMD by smearing the box
grid_spacing = grid_si(2) - grid_si(1);
[~,smear_dir]= max(size(grid_si)); % direction of smearing
dens_profile = smear_phase(dens_profile, psf_DMD, grid_spacing, smear_dir); 
dens_profile(dens_profile < 0) = 0; % make sure that density is positive

end