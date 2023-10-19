clear all
close all
addpath('../classes/')

%gas parameter
temp = 50e-9; %50 nK
n_fourier_cutoff = 15; %the maximum n 
pixel_size = 1e-6; %1 microns 
condensate_length = 100e-6; 
pixnumz = floor(condensate_length/pixel_size); %number of longitudinal pixels
z_grid = linspace(-condensate_length/2, condensate_length/2-pixel_size, pixnumz);
N_samples = 100;

%1. Initiate sampling with Ornstein-Uhlenbeck (OU) process
%Cite references
OU_suite = class_OU_phase_sampling(temp);
rel_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz, N_samples);
com_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz, N_samples);

%com_phase_input_cosineCoeffs = OU_suite.fourier_cosine_coeffs;
%com_phase_input_sineCoeffs = OU_suite.fourier_sine_coeffs;
ext_com_phase_all = zeros(N_samples, pixnumz);
%%%%%Simulating TOF%%%%%%%
count = 0;
for i = 1:N_samples
    interference_suite = class_interference_pattern([rel_phase(i,:);com_phase(i,:)]);

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

    %calculating density ripple
    %the density ripple is not defined in position where mean density is zero
    %(near boundaries), so we introduce a cutoff to the boundary
    boundary_factor = 0.05; %cut-down 5% of pixels from each side of the boundary
    cut1 = boundary_factor*pixnumz;
    cut2 = (1-boundary_factor)*pixnumz;

    ripple = 1 - amp_full(cut1:cut2)./amp_trans(cut1:cut2);

    %Use ripple information to deduce the spectrum of the common phase
    common_suite = class_common_phase_spectrum(ripple, z_grid(cut1:cut2));
    [com_phase_output_cosineCoeffs, com_phase_output_sineCoeffs] = common_suite.extract_com_spectrum(n_fourier_cutoff);
    ext_com_phase = common_suite.extract_com_profile(z_grid);
    ext_com_phase_all(i,:) = ext_com_phase;
    f = common_suite.fidelity_coh(ext_com_phase(cut1:cut2), com_phase(i,cut1:cut2));
    disp(f)
    count = count+1;
    disp(count)

end
