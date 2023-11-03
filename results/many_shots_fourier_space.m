clear all
close all
addpath('../classes/')
addpath('../plotting_func')

%This script is to give an example for single shot common phase and
%relative phase extraction

%Gas parameter
temp_rel = 20e-9;
temp_com = 70e-9;
n_fourier_cutoff = 15 ; %the maximum n 
pixnumz = 100;
t_tof = 15e-3;
N_samples = 500;

flag_error_mitigation = 1;
flag_interaction_broadening = 0;
mean_density_profile = 'BoxPotential'; %Thomas-Fermi mean-density profile

%1. Initiate sampling with Bogoliubov sampling
sampling_suite_rel = class_bogoliubov_phase_sampling(temp_rel);
sampling_suite_com = class_bogoliubov_phase_sampling(temp_com);
rel_phase = sampling_suite_rel.generate_samples(n_fourier_cutoff, pixnumz, N_samples);
com_phase = sampling_suite_com.generate_samples(n_fourier_cutoff, pixnumz, N_samples);

input_cosineCoeffs = sampling_suite_com.fourier_cosine_coeffs;
input_sineCoeffs = sampling_suite_com.fourier_sine_coeffs;

%2. Initializing the output
output_cosineCoeffs = zeros(N_samples, n_fourier_cutoff);
output_sineCoeffs = zeros(N_samples, n_fourier_cutoff);
count =0;
%%%%%3. Simulating TOF%%%%%%%
for i = 1:N_samples
    interference_suite = class_interference_pattern([rel_phase(i,:);com_phase(i,:)], t_tof, mean_density_profile, flag_interaction_broadening);
    x_grid = interference_suite.output_grid_x;
    z_grid = interference_suite.output_grid_z;
    insitu_density = interference_suite.insitu_density';
    condensate_length = interference_suite.condensate_length_Lz;

    rho_tof_full = interference_suite.tof_full_expansion();

    %4. Compute density ripple
    amp_full = trapz(x_grid, rho_tof_full, 2);

    %calculating ripple function
    %the ripple function is not defined in position where mean density is zero
    %(near boundaries), so we introduce a cutoff to the boundary
    boundary_factor = 0.05; %cut-down 5% of pixels from each side of the boundary
    cut1 = floor(boundary_factor*pixnumz);
    cut2 = floor((1-boundary_factor)*pixnumz);

    ripple_func = 1 - amp_full(cut1:cut2)./(2*insitu_density(cut1:cut2));
    %5. Use ripple information to deduce the spectrum of the common phase
    common_suite_raw = class_common_phase_spectrum(ripple_func, z_grid(cut1:cut2), t_tof);
    [raw_output_cosineCoeffs, raw_output_sineCoeffs] = common_suite_raw.extract_com_spectrum(n_fourier_cutoff);
    %6. Error mitigation step
    if flag_error_mitigation == 1
        %extract relative phase
        rel_suite = class_relative_phase_extraction(rho_tof_full, t_tof,flag_interaction_broadening);
        ext_rel_phase = rel_suite.fitting(rel_suite.init_phase_guess());

        %Compute density ripple without common phase
        interference_suite_woc = class_interference_pattern(ext_rel_phase, t_tof, mean_density_profile, flag_interaction_broadening); %woc -> without common phase
        rho_tof_full_woc = interference_suite_woc.tof_full_expansion();

        amp_woc = trapz(x_grid, rho_tof_full_woc, 2);
        ripple_error = 1 - amp_woc(cut1:cut2)./(2*insitu_density(cut1:cut2));
        ripple_func_new = ripple_func - ripple_error;
    else
        ripple_func_new = ripple_func;
    end
    common_suite_mitigated = class_common_phase_spectrum(ripple_func_new, z_grid(cut1:cut2), t_tof);
    [mitigated_cosineCoeffs, mitigated_sineCoeffs] = common_suite_mitigated.extract_com_spectrum(n_fourier_cutoff);
    %if abs(mitigated_cosineCoeffs(1)) < abs(raw_output_cosineCoeffs(1))
        output_cosineCoeffs(i,:) = mitigated_cosineCoeffs;
        output_sineCoeffs(i,:) = mitigated_sineCoeffs;
    %else
     %   output_cosineCoeffs(i,:) = raw_output_sineCoeffs;
     %   output_sineCoeffs(i,:) = raw_output_sineCoeffs;
    %end
   count = count+1;
   disp(count)
end

abs_output_coeffs = zeros(N_samples, n_fourier_cutoff);
abs_input_coeffs = zeros(N_samples, n_fourier_cutoff);
for i = 1:N_samples
    for j = 1:n_fourier_cutoff
        abs_output_coeffs(i,j) = output_cosineCoeffs(i,j)^2+output_sineCoeffs(i,j)^2;
        abs_input_coeffs(i,j) = input_cosineCoeffs(i,j)^2+input_sineCoeffs(i,j)^2;
    end
end

mean_abs_fourier_coeffs_out = zeros(1, n_fourier_cutoff);
mean_abs_fourier_coeffs_in = zeros(1,n_fourier_cutoff);
for i = 1:n_fourier_cutoff
    mean_abs_fourier_coeffs_out(i) = mean(abs_output_coeffs(:,i));
    mean_abs_fourier_coeffs_in(i) = mean(abs_input_coeffs(:,i));
end

save('fourier_stats_70nK_500samples.mat', 'com_phase', 'rel_phase', 'abs_output_coeffs', 'abs_input_coeffs', ...
    'mean_abs_fourier_coeffs_in', 'mean_abs_fourier_coeffs_out','input_sineCoeffs','input_cosineCoeffs', ...
    'output_sineCoeffs', 'output_cosineCoeffs')