function L1 = L1_mat(dens_profile,grid_spacing_si,periodic_BC)

%% some constants

hqsi = 1.054571726e-34; % reduced Planck's constant in SI-units
amuKg  = 1.6605387e-27;         % kg

m = 87*amuKg; % Mass of Rb87 in si units

%% calculate kernel

dens_profile_z_minus_1 = circshift(dens_profile,1,2);
z_times_z_minus_1 = 1./sqrt(dens_profile.*dens_profile_z_minus_1);
z_times_z_minus_1(1) = 0;

diag_part = 2*diag(1./dens_profile);
off_diag_part = -(diag(z_times_z_minus_1(2:end),1) + diag(z_times_z_minus_1(2:end),-1));

%------ BC
if periodic_BC
    off_diag_part(1,end) = - 1/sqrt(dens_profile(1)*dens_profile(end));
    off_diag_part(end,1) = off_diag_part(1,end);
end
%-------

L1 = 1/grid_spacing_si * hqsi^2/(2*m) * (diag_part + off_diag_part);

end