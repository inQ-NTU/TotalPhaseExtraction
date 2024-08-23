close all; clear all;
addpath('../classes/')
addpath('../plotting_func/')

%set parameter
T_com = 40e-9;
T_rel = 40e-9;
mean_density = 75e6;
condensate_length = 100e-6;
buffer_length = [0, 10e-6]; % [buffer z, buffer x]
pixnumz = (condensate_length*1e6)+1;
t_tof = 11e-3;
%num_samples = 10000;
max_fourier = 40;
%z_grid = linspace(-condensate_length/2, condensate_length/2, pixnumz);

%begin sampling
com_sampling_suite = class_bogoliubov_sampling(T_com, mean_density, condensate_length);
rel_sampling_suite = class_bogoliubov_sampling(T_rel, mean_density, condensate_length);
[com_phase_sample, com_density_fluct] = com_sampling_suite.generate_fluct_samples(max_fourier, pixnumz);
[rel_phase_sample, rel_density_fluct] = rel_sampling_suite.generate_fluct_samples(max_fourier, pixnumz);

phi_1 = (com_phase_sample + rel_phase_sample)/2;
phi_2 = (com_phase_sample - rel_phase_sample)/2;

density_fluct_1 = (com_density_fluct + rel_density_fluct)/2;
density_fluct_2 = (com_density_fluct - rel_density_fluct)/2;

insitu_density_1 = mean_density+density_fluct_1;
insitu_density_2 = mean_density+density_fluct_2;

interference_suite = class_interference_pattern_new([phi_1; phi_2], [insitu_density_1; insitu_density_2], t_tof);
rho_tof_trans = interference_suite.tof_transverse_expansion();
rho_tof_full = interference_suite.tof_full_expansion();

x_grid = interference_suite.output_grid_x;
z_grid = interference_suite.output_grid_z;

figure
imagesc(x_grid, z_grid, rho_tof_trans)
colorbar
colormap(gge_colormap)

figure
imagesc(x_grid, z_grid, rho_tof_full)
colorbar
colormap(gge_colormap)

rel_suite_trans = class_relative_phase_extraction(rho_tof_trans, t_tof, x_grid);
rel_suite_full = class_relative_phase_extraction(rho_tof_full, t_tof, x_grid);

ext_rel_trans = rel_suite_trans.fitting(rel_suite_trans.init_phase_guess());
ext_rel_full = rel_suite_full.fitting(rel_suite_full.init_phase_guess());

figure
plot(z_grid, rel_phase_sample)
hold on
plot(z_grid, ext_rel_trans, 'o')
plot(z_grid, ext_rel_full, 'x-')



