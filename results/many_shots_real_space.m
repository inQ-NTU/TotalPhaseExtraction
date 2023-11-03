clear all
close all
addpath('../classes/')
addpath('../plotting_func')

%This script is to give an example for single shot common phase and
%relative phase extraction

%Gas parameter
temp_com = 50e-9;
temp_rel = 20e-9;
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

%2. Initializing the output
raw_ext_com_phase = zeros(N_samples, pixnumz);
mitigated_ext_com_phase = zeros(N_samples, pixnumz);
raw_fidelity = zeros(1,N_samples);
mitigated_fidelity = zeros(1,N_samples);

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
    raw_output = common_suite_raw.extract_com_profile(z_grid);
    raw_ext_com_phase(i,:) = raw_output;
    raw_fid = common_suite_raw.fidelity_coh(com_phase(i,:), raw_output);
    disp(raw_fid)
    raw_fidelity(i) =raw_fid;
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
    [output_cosineCoeffs, output_sineCoeffs] = common_suite_mitigated.extract_com_spectrum(n_fourier_cutoff);
    %if abs(output_cosineCoeffs(1)) < abs(raw_output_cosineCoeffs(1))
        mitigated_output = common_suite_mitigated.extract_com_profile(z_grid);
        mitigated_ext_com_phase(i,:) = mitigated_output;
        fid = common_suite_mitigated.fidelity_coh(com_phase(i,:), mitigated_output);
        disp(fid)
        mitigated_fidelity(i) = fid;
    %else
     %   mitigated_ext_com_phase(i,:) = raw_output;
     %   mitigated_fidelity(i) = raw_fid;
    %end
   count = count+1;
   disp(count)
end


input_corr = class_1d_correlation(com_phase);
cov_in = input_corr.covariance_matrix();

final_corr = class_1d_correlation(mitigated_ext_com_phase);
cov_out = final_corr.covariance_matrix();

save('real_space_reconstruction_50nk_500samples.mat', 'com_phase', 'rel_phase', 'raw_ext_com_phase', 'raw_fidelity', ...
    'mitigated_ext_com_phase','mitigated_fidelity', 'cov_in', 'cov_out')


if 0
load('many_shots_real_space.mat')
condensate_length = 100; %100 microns
z_grid = linspace(-condensate_length/2,condensate_length/2,pixnumz);

figure
g = tight_subplot(1,2,[.12 .12],[.1 .05],[.1 .05]);

axes(g(1))
histogram(raw_fidelity,'Normalization','probability')
hold on
histogram(mitigated_fidelity,'Normalization','probability')
xlim([0.5,1])
xticks([0.5,0.6,0.7,0.8,0.9,1])
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.25,0.8]);
ylim([0,0.21])
yticks([0,0.1,0.2])

axes(g(2))
plot(z_grid, cov_in(50,:),'Color','Black', 'LineStyle', '-.')
hold on
plot(z_grid, cov_out(50,:))
end


