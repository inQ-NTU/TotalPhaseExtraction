function [atoms,width] = func_create_atomic_density...
    (psi_L,psi_R,TOF_si,grid_si,fringeSpacing,omega_trans,transversal_flag)

% transversal_flag.... 0 for VAndor images
%                      1 for TAndor

%% constants
hbar     = 6.626e-34/(2*pi);      % Js
amuKg    = 1.6605387e-27;         % kg
Aatom = 87;                 % Rb atoms
massAtom = amuKg*Aatom;

%% calculate atomic densities

% free expansion in 1D
psi_L_TOF = func_free_TOF_propagation_1D(psi_L',grid_si,TOF_si); 
psi_R_TOF = func_free_TOF_propagation_1D(psi_R',grid_si,TOF_si); 

% free expansion transversal
width0=sqrt(hbar/(massAtom*omega_trans)); %in trap width = harmonic oscillator ground state width
width=width0*sqrt(1+(TOF_si*omega_trans)^2); %after free expansion

halfd = pi * massAtom * width0^2 * width^2 /(hbar * TOF_si * fringeSpacing); % half the distance between the wells

%calculate 1D densities
psi_L_TOF_density = conj(psi_L_TOF).*psi_L_TOF; %left
psi_R_TOF_density = conj(psi_R_TOF).*psi_R_TOF; %right
coh_TOF=psi_L_TOF.*conj(psi_R_TOF); %interference term

if transversal_flag %for artificial TAndor images    
    atoms1 = 1/(sqrt(pi)*width) * exp(-(grid_si').^2/width^2) * psi_L_TOF_density; %left    
    atoms2 = 1/(sqrt(pi)*width) * exp(-(grid_si').^2/width^2) * psi_R_TOF_density; %right
    atoms3 = zeros(size(atoms2)); %interference term
else %for artificial VAndor images
    atoms1 = 1/(sqrt(pi)*width) * exp(-(grid_si'+halfd).^2/width^2) * psi_L_TOF_density;
    atoms2 = 1/(sqrt(pi)*width) * exp(-(grid_si'-halfd).^2/width^2) * psi_R_TOF_density;
    atoms3 = 2*real( 1/(sqrt(pi)*width) * (exp(-(grid_si'.^2 + halfd^2)/width^2) .* exp(1i*2*pi*grid_si'/fringeSpacing)) * coh_TOF); %interference term
end
    
atoms = atoms1 + atoms2 + atoms3;

end