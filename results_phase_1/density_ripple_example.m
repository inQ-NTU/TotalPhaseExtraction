clear all
close all
addpath('../classes/')
addpath('../plotting_func')

%This script is to give an example for single shot common phase and
%relative phase extraction

%Gas parameter
temp = 40e-9; 
n_fourier_cutoff = 20; %the maximum n 
pixnumz = 100;
t_tof = 15e-3;

flag_error_mitigation = 1;
flag_interaction_broadening = 0;
mean_density_profile = 'BoxPotential'; %Thomas-Fermi mean-density profile

%if 0
%1. Initiate sampling with Bogoliubov sampling
sampling_suite = class_bogoliubov_phase_sampling(temp);
rel_phase = sampling_suite.generate_samples(n_fourier_cutoff, pixnumz);
com_phase = sampling_suite.generate_samples(n_fourier_cutoff, pixnumz);

%2. Fourier analysis of the input common phase
input_cosineCoeffs = sampling_suite.fourier_cosine_coeffs;
input_sineCoeffs = sampling_suite.fourier_sine_coeffs;
%end

%load('input_single_shot_phase.mat')
%%%%%3. Simulating TOF%%%%%%%
interference_suite = class_interference_pattern([rel_phase;com_phase], t_tof, mean_density_profile, flag_interaction_broadening);
interference_suite_woc = class_interference_pattern(rel_phase, t_tof, mean_density_profile, flag_interaction_broadening);

x_grid = interference_suite.output_grid_x;
z_grid = interference_suite.output_grid_z;
insitu_density = interference_suite.insitu_density';
condensate_length = interference_suite.condensate_length_Lz;

rho_tof_trans = interference_suite.tof_transversal_expansion();
rho_tof_full = interference_suite.tof_full_expansion();
rho_tof_full_woc = interference_suite_woc.tof_full_expansion();

amp_full = trapz(x_grid, rho_tof_full, 2);
amp_trans = trapz(x_grid, rho_tof_trans, 2);
amp_full_woc = trapz(x_grid, rho_tof_full_woc, 2);

x_grid = linspace(-60,60,120);
z_grid = linspace(-50,50,100);

g = tight_subplot(2,1,[.08 .1],[.2 .2],[.2 .2]);
axes(g(1))
imagesc(z_grid, x_grid, rho_tof_full'.*1e-12)
colormap(gge_colormap)
clim([0,5])
cb = colorbar('Location','northoutside');
cb.Position = cb.Position + [0,0.1,0,0];
ylabel(cb, '$\rho\; (\rm\mu m^{-2})$', 'Interpreter','latex','FontSize',20)
xticks([])
yticks([-60,-20,20,60])
ylabel('$x \;(\rm \mu m)$', 'Interpreter','latex')


axes(g(2))
plot(z_grid, amp_full.*1e-6, 'Color','red')
hold on
plot(z_grid, amp_trans.*1e-6, 'Color', 'black','LineStyle', '-.')
plot(z_grid, amp_full_woc.*1e-6, 'Color','red','LineStyle','--')
ylabel('$n\; (\rm\mu m^{-1})$', 'Interpreter','latex')
xlabel('$z \;(\rm \mu m)$', 'Interpreter','latex')

set(g, 'FontName', 'Times', 'FontSize', 20)