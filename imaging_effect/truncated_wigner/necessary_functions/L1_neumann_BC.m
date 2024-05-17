function L1 = L1_neumann_BC(dens_profile,grid_spacing_si)

%% some constants

hqsi = 1.054571726e-34; % reduced Planck's constant in SI-units
amuKg  = 1.6605387e-27;         % kg

m = 87*amuKg; % Mass of Rb87 in si units

%% calculate kernel

dens_profile_z_minus_1 = circshift(dens_profile,1,2);
z_times_z_minus_1 = 1./sqrt(dens_profile.*dens_profile_z_minus_1);
z_times_z_minus_1(1) = 0;

diag_part = 2*diag(1./dens_profile);

% neumann BC, derivative zero at the boundary
diag_part(1,1) = diag_part(1,1)/2;
diag_part(end,end) = diag_part(end,end)/2;

off_diag_part = -(diag(z_times_z_minus_1(2:end),1) + diag(z_times_z_minus_1(2:end),-1));

L1 = 1/grid_spacing_si * hqsi^2/(2*m) * (diag_part + off_diag_part);

end