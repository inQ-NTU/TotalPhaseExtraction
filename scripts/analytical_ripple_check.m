clear all
close all
addpath('../classes/')
addpath('../plotting_func')

%This script is to give an example for single shot common phase and
%relative phase extraction

%gas parameter
temp = 50e-9; %50 nK
n_fourier_cutoff = 12; %the maximum n 
pixel_size = 1e-6; %1 microns 
condensate_length = 100e-6; 
pixnumz = floor(condensate_length/pixel_size); %number of longitudinal pixels
z_grid = linspace(-condensate_length/2, condensate_length/2-pixel_size, pixnumz);
x_grid = z_grid;

flag_interaction_broadening = 1; %turn-on/off interaction broadening
mean_density_profile = 'BoxPotential'; %Thomas-Fermi mean-density profile
t_tof = 15e-3;

%1. Initiate sampling with Ornstein-Uhlenbeck (OU) process
%Cite references
%if 0
OU_suite = class_OU_phase_sampling(temp);
%rel_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz);
rel_phase = zeros(1,pixnumz);
%com_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz);
a = 1;
n = 3;
com_phase = acos(n*2*pi*z_grid);
eps_t = sqrt(OU_suite.hbar*t_tof/OU_suite.m)/condensate_length;
analytical_ripple = -(2*pi*n)^2*a*((eps_t^2)/2)*cos(n*2*pi*z_grid/condensate_length);
%disp(eps_t)
if 0
input_cosineCoeffs = OU_suite.fourier_cosine_coeffs;
input_sineCoeffs = OU_suite.fourier_sine_coeffs;

analytical_ripple = zeros(1, pixnumz);
for i = 1:pixnumz
    val = 0;
    for n = 1:n_fourier_cutoff
        val = val - (2*pi*n)^2*(input_cosineCoeffs(n)*cos(2*pi*n*z_grid(i)/condensate_length)+...
            input_sineCoeffs(n)*sin(2*pi*n*z_grid(i)/condensate_length));
    end
    analytical_ripple(i) = val;
end
analytical_ripple = analytical_ripple*eps_t^2/2;
end
%end
%load('single_shot_phase_profile.mat')
%%%%%Simulating TOF%%%%%%%
interference_suite = class_interference_pattern([rel_phase;com_phase],t_tof, 'BoxPotential', flag_interaction_broadening);

rho_tof_trans = interference_suite.tof_transversal_expansion();
rho_tof_full = interference_suite.tof_full_expansion();

%Normalize to 10,000 atoms
rho_tof_trans = interference_suite.normalize(rho_tof_trans, 10^4);
rho_tof_full = interference_suite.normalize(rho_tof_full, 10^4);

%Compute density ripple
amp_trans = trapz(x_grid, rho_tof_trans, 2);
amp_full = trapz(x_grid, rho_tof_full, 2);

%calculating ripple function
%the ripple function is not defined in position where mean density is zero
%(near boundaries), so we introduce a cutoff to the boundary
boundary_factor = 0.1; %cut-down 10% of pixels from each side of the boundary
cut1 = boundary_factor*pixnumz;
cut2 = (1-boundary_factor)*pixnumz;

ripple_func = 1 - amp_full(cut1:cut2)./amp_trans(cut1:cut2);
common_suite = class_common_phase_spectrum(ripple_func, z_grid(cut1:cut2), t_tof);
[output_cosineCoeffs, output_sineCoeffs] = common_suite.extract_com_spectrum(n_fourier_cutoff);
ext_com_phase = common_suite.extract_com_profile(z_grid);

ft = fft(analytical_ripple)*0.01;

plot(z_grid(cut1:cut2),ripple_func)
hold on
plot(z_grid, analytical_ripple)
