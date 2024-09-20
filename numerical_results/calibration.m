clear all; close all
addpath('../classes/')

%Define gas and extraction parameters as well as simulation parameters
temp_com = 50e-9; %temperature of the common phase (nK)
temp_rel = 30e-9; %temperature of the relative phase (nK)
t_tof_1 = 11e-3; %short expansion time
t_tof_2 = 16e-3; %long expansion_time
mean_density = 75e6; %Background density of each 1D gas (atoms/m)

%By default, the gas lies on the interval[-50 microns, 50 microns] insitu 
bulk_start = -40e-6; %Only extract common phase starting from this position (ignore edges)
bulk_end = 40e-6; %Only extract common phase until this position (ignore edges)

%Simulation parameter
n_fourier_cutoff = 40; %maximum mode num in sampling the common phase
ext_cutoff = 40; %maximum mode in the extraction 

pixnumz = 201; %number of longitudinal points
num_samples = 10; %number of phase samples
scaling_factor = sqrt(2);

%Flags
flag_density_fluct = 1;

%Sample in situ fluctuations with Bogoliubov sampling
sampling_suite_rel = class_bogoliubov_sampling(temp_rel, mean_density, scaling_factor);
sampling_suite_com = class_bogoliubov_sampling(temp_com, mean_density, scaling_factor);

%Generate fluctuation samples
%[rel_phase, rel_density_fluct] = sampling_suite_rel.generate_fluct_samples(n_fourier_cutoff, pixnumz, num_samples);
%[com_phase, com_density_fluct] = sampling_suite_com.generate_fluct_samples(n_fourier_cutoff, pixnumz, num_samples);
rel_phase = zeros(num_samples, pixnumz);
com_phase = zeros(num_samples,pixnumz);
com_density_fluct = zeros(num_samples,pixnumz);
rel_density_fluct = zeros(num_samples,pixnumz);
%Calculate input fluctuation for each gas
phase_1 = (com_phase + rel_phase)/2;
phase_2 = (com_phase - rel_phase)/2;

density_fluct_1 = (com_density_fluct + rel_density_fluct)/2;
density_fluct_2 = (com_density_fluct - rel_density_fluct)/2;

count = 0;
%%%%%Simulating TOF and collecting density ripple data%%%%%%%
for i = 1:num_samples
    if flag_density_fluct
        insitu_density_1 = mean_density*ones(1,pixnumz) + density_fluct_1(i,:); %adding density fluctuation into mean density
        insitu_density_2 = mean_density*ones(1,pixnumz) + density_fluct_2(i,:); %adding density fluctuation into mean density
    else
        insitu_density_1 = mean_density*ones(1,pixnumz);
        insitu_density_2 = mean_density*ones(1,pixnumz);
    end
    interference_suite_t1 = class_interference_pattern([phase_1(i,:); phase_2(i,:)], [insitu_density_1; insitu_density_2],...
        t_tof_1);
    interference_suite_t2 = class_interference_pattern([phase_1(i,:); phase_2(i,:)], [insitu_density_1; insitu_density_2],...
        t_tof_2);

    grid_dens = interference_suite_t1.output_grid_z; %grid in longitudinal direction
    x_grid = interference_suite_t1.output_grid_x; 

    rho_tof_full_t1 = interference_suite_t1.tof_full_expansion();
    rho_tof_full_t2 = interference_suite_t2.tof_full_expansion();

    amp_full_t1 = trapz(x_grid, rho_tof_full_t1, 2);
    amp_full_t2 = trapz(x_grid, rho_tof_full_t2, 2);

    density_ripple_t1(i,:) = amp_full_t1;
    density_ripple_t2(i,:) = amp_full_t2;

    %Ignoring the edges
    [val, idx_1] = min(abs(grid_dens - bulk_start)); %finding start index of the bulk
    [val, idx_2] = min(abs(grid_dens - bulk_end)); % finding end index of the bulk
    cut_grid_dens = grid_dens(idx_1:idx_2); %truncated z grid

    cut_amp_full_t1 = amp_full_t1(idx_1:idx_2);
    cut_amp_full_t2 = amp_full_t2(idx_1:idx_2);
    cut_density_ripple_t1(i,:) = cut_amp_full_t1;
    cut_density_ripple_t2(i,:) = cut_amp_full_t2;
    
    %Extract common phase
    common_suite_t1 = class_common_phase_extraction(cut_amp_full_t1, 2*mean_density*ones(1,length(cut_grid_dens)), ...
        cut_grid_dens, t_tof_1);
    common_suite_t2 = class_common_phase_extraction(cut_amp_full_t2, 2*mean_density*ones(1,length(cut_grid_dens)), ...
        cut_grid_dens, t_tof_2);
    [cosineCoeffs_t1, sineCoeffs_t1] = common_suite_t1.extract_com_spectrum(ext_cutoff);
    [cosineCoeffs_t2, sineCoeffs_t2] = common_suite_t2.extract_com_spectrum(ext_cutoff);
    out_cosineCoeffs_t1(i,:) = cosineCoeffs_t1;
    out_cosineCoeffs_t2(i,:) = cosineCoeffs_t2;
    out_sineCoeffs_t1(i,:) = sineCoeffs_t1;
    out_sineCoeffs_t2(i,:) = sineCoeffs_t2;
    out_com_phase_t1(i,:) = common_suite_t1.reconstruct_com_phase(cosineCoeffs_t1, sineCoeffs_t1);
    out_com_phase_t2(i,:) = common_suite_t2.reconstruct_com_phase(cosineCoeffs_t2, sineCoeffs_t2);

    count = count+1
end

spectrum_t1 = out_cosineCoeffs_t1.^2+out_sineCoeffs_t1.^2;
spectrum_t2 = out_cosineCoeffs_t2.^2+out_sineCoeffs_t2.^2;
mean_spectrum_t1 = mean(spectrum_t1, 1);
mean_spectrum_t2 = mean(spectrum_t2, 1);

%plot(cut_grid_dens, out_com_phase_t1(1,:))
%hold on
%plot(cut_grid_dens, out_com_phase_t2(1,:))

%save('scan_11ms_16ms_50nk_LE.mat', 't_tof_1', 't_tof_2', 'temp_com', 'temp_rel', 'mean_density', 'n_fourier_cutoff', 'ext_cutoff', 'grid_dens','cut_grid_dens', 'out_com_phase_t1', ...
%   'out_com_phase_t2','com_phase', 'rel_phase', 'com_density_fluct','rel_density_fluct', 'flag_density_fluct', 'num_samples', 'pixnumz', 'density_ripple_t1', 'cut_density_ripple_t1',...
%   'density_ripple_t2', 'cut_density_ripple_t2', 'idx_1', 'idx_2', 'out_sineCoeffs_t2', 'out_sineCoeffs_t1', 'out_cosineCoeffs_t2', 'out_cosineCoeffs_t1')