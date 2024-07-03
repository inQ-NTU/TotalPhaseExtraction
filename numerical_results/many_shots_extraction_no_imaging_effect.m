clear all; close all
addpath('../classes/')

%Define gas and extraction parameters as well as simulation parameters
temp_com = 50e-9; %temperature of the common phase (nK)
temp_rel = 30e-9; %temperature of the relative phase (nK)
t_tof = 11e-3; %expansion time
mean_density = 75e6; %Background density of each 1D gas (atoms/microns)

%By default, the gas lies on the interval[-50 microns, 50 microns] insitu 
bulk_start = -40e-6; %Only extract common phase starting from this position (ignore edges)
bulk_end = 40e-6; %Only extract common phase until this position (ignore edges)

%Simulation parameter
n_fourier_cutoff = 40; %maximum mode num in sampling the common phase
ext_cutoff = 40; %maximum mode in the extraction 
pixnumz = 201; %number of longitudinal points
num_samples = 1; %number of phase samples

%Flags
flag_density_fluct = 1;

%Initiate sampling with Bogoliubov sampling
%load('input_comparison_imaging.mat')
sampling_suite_rel = class_bogoliubov_sampling(temp_rel, mean_density);
sampling_suite_com = class_bogoliubov_sampling(temp_com, mean_density);
sampling_suite_density = class_bogoliubov_sampling(temp_com, mean_density);

%Generate fluctuation samples
rel_phase = sampling_suite_rel.generate_phase_samples(n_fourier_cutoff, pixnumz, num_samples);
com_phase = sampling_suite_com.generate_phase_samples(n_fourier_cutoff, pixnumz, num_samples);
density_fluct = sampling_suite_com.generate_density_fluct_samples(n_fourier_cutoff, pixnumz, num_samples);

%for calibration
%rel_phase = zeros(num_samples, pixnumz);
%com_phase = zeros(num_samples, pixnumz);
%density_fluct = zeros(num_samples, pixnumz);

count = 0;
%%%%%Simulating TOF and collecting density ripple data%%%%%%%
for i = 1:num_samples
    if flag_density_fluct
        insitu_density = mean_density + density_fluct(i,:); %adding density fluctuation into mean density
    else
        insitu_density = mean_density;
    end
    interference_suite = class_interference_pattern([rel_phase(i,:);com_phase(i,:)], t_tof, insitu_density);
    condensate_length = interference_suite.condensate_length_Lz;
    grid_dens = interference_suite.output_grid_z; %grid in longitudinal direction
    x_grid = interference_suite.output_grid_x; 

    rho_tof_full = interference_suite.tof_full_expansion();

    amp_full = trapz(x_grid, rho_tof_full, 2);
    density_ripple(i,:) = amp_full;

    %Ignoring the edges
    [val, idx_1] = min(abs(grid_dens - bulk_start)); %finding start index of the bulk
    [val, idx_2] = min(abs(grid_dens - bulk_end)); % finding end index of the bulk
    cut_grid_dens = grid_dens(idx_1:idx_2); %truncated z grid

    cut_amp_full =   amp_full(idx_1:idx_2);
    cut_density_ripple(i,:) = cut_amp_full;
    cut_condensate_length = abs(bulk_end-bulk_start);
    
    %Extract common phase
    common_suite = class_common_phase_extraction(cut_amp_full, 2*mean_density*ones(1,length(cut_grid_dens)), ...
        cut_grid_dens,cut_condensate_length, t_tof);
    [cosineCoeffs_u, sineCoeffs_u] = common_suite.extract_com_spectrum_uniform(ext_cutoff);
    com_phase_u = common_suite.reconstruct_com_phase(cosineCoeffs_u, sineCoeffs_u);
    output_com_phase(i,:) = com_phase_u;

    count = count+1
end
%save('scan_11ms_calibrate_LE_RP_DF.mat', 't_tof', 'temp_com', 'temp_rel', 'mean_density', 'n_fourier_cutoff', 'ext_cutoff', 'grid_dens','cut_grid_dens', 'output_com_phase', ...
%   'com_phase', 'rel_phase', 'density_fluct', 'flag_density_fluct', 'num_samples', 'pixnumz', 'density_ripple', 'cut_density_ripple')