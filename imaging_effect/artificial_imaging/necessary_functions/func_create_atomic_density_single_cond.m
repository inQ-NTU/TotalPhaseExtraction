function [atoms,width] = func_create_atomic_density_single_cond(psi,TOF_si,grid_si,omega_trans)

% transversal_flag.... 0 for VAndor images
%                      1 for TAndor

%% constants
hbar     = 6.626e-34/(2*pi);      % Js
amuKg    = 1.6605387e-27;         % kg
Aatom = 87;                 % Rb atoms
massAtom = amuKg*Aatom;

%% calculate atomic densities

% free expansion in 1D
psi_TOF = func_free_TOF_propagation_1D(psi',grid_si,TOF_si); 

% free expansion transversal
width0=sqrt(hbar/(massAtom*omega_trans)); %in trap width = harmonic oscillator ground state width
width=width0*sqrt(1+(TOF_si*omega_trans)^2); %after free expansion

% calculate 1D densities
psi_TOF_density = conj(psi_TOF).*psi_TOF; %left

% 2D atomic density
atoms = 1/(sqrt(pi)*width) * exp(-(grid_si').^2/width^2) * psi_TOF_density; %left

end