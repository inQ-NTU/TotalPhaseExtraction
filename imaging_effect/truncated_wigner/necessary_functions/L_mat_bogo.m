function L_mat = L_mat_bogo(dens_profile,g0_si_arr,J_Hz,grid_spacing_si,periodic_BC)

%% some constants

hqsi = 1.054571726e-34; % reduced Planck's constant in SI-units
amuKg  = 1.6605387e-27;         % kg

m = 87*amuKg; % Mass of Rb87 in si units

%% calculate kernel

if periodic_BC
    L1 = L1_mat(dens_profile,grid_spacing_si,periodic_BC);
else
    L1 = L1_neumann_BC(dens_profile,grid_spacing_si);
end

L2 = L2_mat(dens_profile,grid_spacing_si,periodic_BC);
L3 = 2 * grid_spacing_si * diag(g0_si_arr);
Lt = 2 * grid_spacing_si * hqsi * diag(2*pi*J_Hz./dens_profile);

L_mat = L1 + L2 + L3 + Lt;

%just for testing
% L_mat = L1 + L3 + Lt;

end