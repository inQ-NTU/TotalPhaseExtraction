clear all; close all
addpath('../classes/')
addpath('../plotting_func')

%Define gas and extraction parameters as well as simulation parameters
temp_com = 30e-9; %temperature of the common phase (nK)
temp_rel = 10e-9; %temperature of the relative phase (nK)
t_tof = 20e-3; %expansion time
mean_density = 75e6; %Background density of each 1D gas (atoms/microns)

%By default, the gas lies on the interval[-50 microns, 50 microns] insitu 
buffer_length = 10e-6; %10 microns buffer for long expansion time (>20 ms)
bulk_start = -40e-6; %Only extract common phase starting from this position (ignore edges)
bulk_end = 40e-6; %Only extract common phase until this position (ignore edges)

%Simulation parameter
n_fourier_cutoff = 40; %maximum mode num in sampling the common phase
ext_cutoff = 40; %maximum mode in the extraction 
pixnumz = 200; %number of longitudinal points
num_samples = 1000; %number of phase samples

%Flags
flag_density_fluct = 1;
flag_buffer = 0;

%Initiate sampling with Bogoliubov sampling
sampling_suite_rel = class_bogoliubov_sampling(temp_rel, mean_density);
sampling_suite_com = class_bogoliubov_sampling(temp_com, mean_density);

%Generate fluctuation samples
rel_phase = sampling_suite_rel.generate_phase_samples(n_fourier_cutoff, pixnumz, num_samples);
com_phase = sampling_suite_com.generate_phase_samples(n_fourier_cutoff, pixnumz, num_samples);
density_fluct = sampling_suite_com.generate_density_fluct_samples(n_fourier_cutoff, pixnumz, num_samples);

%%%%%3. Simulating TOF%%%%%%%
count = 0;
for i = 1:num_samples
    if flag_density_fluct
        insitu_density = mean_density + density_fluct(i,:); %adding density fluctuation into mean density
    else
        insitu_density = mean_density;
    end

    if flag_buffer
        interference_suite = class_interference_pattern([rel_phase(i,:);com_phase(i,:)], t_tof, insitu_density, buffer_length);
    else
        interference_suite = class_interference_pattern([rel_phase(i,:);com_phase(i,:)], t_tof, insitu_density);
    end
    
    x_grid = interference_suite.output_grid_x; %grid in transverse direction (possibly buffered)
    buffered_z_grid = interference_suite.output_grid_z; %grid in longitudinal direction (possibly buffered)

    %Unbuffered longitudinal grid
    condensate_length = interference_suite.condensate_length_Lz;
    z_grid = linspace(-condensate_length/2, condensate_length/2, pixnumz);

    [val, idx_1] = min(abs(z_grid - bulk_start)); %finding start index of the bulk
    [val, idx_2] = min(abs(z_grid - bulk_end)); % finding end index of the bulk
    cut_z_grid = z_grid(idx_1:idx_2); %truncated z grid
    
    insitu_mean_density = 2*mean_density*ones(1,length(cut_z_grid));

    rho_tof_full = interference_suite.tof_full_expansion();
    
    amp_full = trapz(x_grid, rho_tof_full, 2);
    [val, idx_1_buffered] = min(abs(buffered_z_grid - bulk_start));
    [val, idx_2_buffered] = min(abs(buffered_z_grid - bulk_end));
    amp_full = amp_full(idx_1_buffered:idx_2_buffered);

    %6. Extract common phase
    common_suite = class_common_phase_extraction(amp_full', insitu_mean_density, cut_z_grid, condensate_length,t_tof);
    [cosineCoeffs_u, sineCoeffs_u] = common_suite.extract_com_spectrum_uniform(ext_cutoff);
    com_phase_u = common_suite.reconstruct_com_phase(cosineCoeffs_u, sineCoeffs_u);
    output_com_phase(i,:) = com_phase_u;
    output_fidelity(i) = common_suite.fidelity_coh(com_phase_u, com_phase(i,idx_1:idx_2))
end

sxave('scan_20ms_30nK', 't_tof', 'temp_com', 'temp_rel', 'mean_density', 'n_fourier_cutoff', 'ext_cutoff', 'cut_z_grid', 'z_grid', 'output_com_phase', 'output_fidelity', ...
    'com_phase', 'rel_phase', 'density_fluct', 'flag_buffer', 'flag_density_fluct', 'num_samples', 'pixnumz')