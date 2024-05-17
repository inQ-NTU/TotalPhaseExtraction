function density_fluct_profiles = bogo_density_fluct_shots(T_si,dens_profile,g0_si_arr,J_Hz_arr,grid_spacing_si,periodic_BC,no_real)

kbsi = 1.3806488e-23; %Boltzmann's constant in SI-units

beta = 1/(kbsi*T_si);

L_mat = L_mat_bogo(dens_profile,g0_si_arr,J_Hz_arr,grid_spacing_si,periodic_BC);
kernel = beta*L_mat;
density_fluct_profiles = shots_from_gaussian(kernel,no_real);

end