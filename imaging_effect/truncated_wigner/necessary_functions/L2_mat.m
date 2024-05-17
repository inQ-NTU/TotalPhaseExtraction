function L2 = L2_mat(dens_profile,grid_spacing_si,periodic_BC)

%% some constants

hqsi = 1.054571726e-34; % reduced Planck's constant in SI-units
amuKg  = 1.6605387e-27;         % kg

m = 87*amuKg; % Mass of Rb87 in si units

%% calculate kernel

dens_profile_z_plus_1 = circshift(dens_profile,-1,2);
dens_profile_z_minus_1 = circshift(dens_profile,1,2);

second_der = 1/grid_spacing_si^2 * (sqrt(dens_profile_z_plus_1) - 2*sqrt(dens_profile) + sqrt(dens_profile_z_minus_1));

if ~periodic_BC
    second_der(1) = second_der(2);
    second_der(end) = second_der(end-1);
    
    %just for testing
%     second_der(1) = 1/grid_spacing_si^2 * (sqrt(dens_profile_z_plus_1(1)) - sqrt(dens_profile(1)));
%     second_der(end) = 1/grid_spacing_si^2 * (sqrt(dens_profile_z_minus_1(end)) - sqrt(dens_profile(end)));
end

L2 = grid_spacing_si * hqsi^2/(2*m) * diag(second_der./(dens_profile.^(3/2)));

end