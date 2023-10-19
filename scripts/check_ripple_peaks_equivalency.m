clear all
close all
addpath('../classes/')
addpath('../plotting_func')

%This script is to give an example for single shot common phase and
%relative phase extraction


%gas parameter
temp = 50e-9; %50 nK
n_fourier_cutoff = 15; %the maximum n 
pixel_size = 1e-6; %1 microns 
condensate_length = 100e-6; 
pixnumz = floor(condensate_length/pixel_size); %number of longitudinal pixels
z_grid = linspace(-condensate_length/2, condensate_length/2-pixel_size, pixnumz);
x_grid = z_grid;

%1. Initiate sampling with Ornstein-Uhlenbeck (OU) process
%Cite references
OU_suite = class_OU_phase_sampling(temp);
rel_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz);
com_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz);

com_phase_input_cosineCoeffs = OU_suite.fourier_cosine_coeffs;
com_phase_input_sineCoeffs = OU_suite.fourier_sine_coeffs;

%%%%%Simulating TOF%%%%%%%
interference_suite = class_interference_pattern([rel_phase;com_phase]);

rho_tof_trans = interference_suite.tof_transversal_expansion();
rho_tof_full = interference_suite.tof_full_expansion();

%initiate relative phase extraction
rel_ext_suite_full = class_relative_phase_extraction(rho_tof_full);
rel_ext_suite_trans = class_relative_phase_extraction(rho_tof_trans);

%extracting relative phase
ext_rel_phase_trans = rel_ext_suite_trans.fitting(rel_ext_suite_trans.init_phase_guess());
ext_rel_phase_full = rel_ext_suite_full.fitting(rel_ext_suite_full.init_phase_guess());

%extracting interference peak
amp_trans = rel_ext_suite_trans.normalization_amplitudes;
amp_full = rel_ext_suite_full.normalization_amplitudes;

%Computing density ripple
ripple_trans = trapz(x_grid, rho_tof_trans, 2);
ripple_full = trapz(x_grid, rho_tof_full, 2);

%plotting
plot(ripple_trans)
hold on
plot(ripple_full)

figure
plot(amp_trans)
hold on
plot(amp_full)