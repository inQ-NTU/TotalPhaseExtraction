clear all; close all
addpath('../classes/')
addpath('../plotting_func')
addpath('../imaging_effect')
addpath('../imaging_effect/utility')
addpath('../imaging_effect/image_analysis')
addpath('../imaging_effect/artificial_imaging')
addpath('../imaging_effect/artificial_imaging/necessary_functions')
addpath('main_Data/')

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
num_samples = 1000; %number of phase samples

%Flags
flag_density_fluct = 1;
scaling_factor = sqrt(2);

%Sample in situ fluctuations with Bogoliubov sampling
%sampling_suite_rel = class_bogoliubov_sampling(temp_rel, mean_density, scaling_factor);
%sampling_suite_com = class_bogoliubov_sampling(temp_com, mean_density, scaling_factor);

%Generate fluctuation samples
%[rel_phase, rel_density_fluct] = sampling_suite_rel.generate_fluct_samples(n_fourier_cutoff, pixnumz, num_samples);
%[com_phase, com_density_fluct] = sampling_suite_com.generate_fluct_samples(n_fourier_cutoff, pixnumz, num_samples);

%[in_cosineCoeffs, in_sineCoeffs] = sampling_suite_com.compute_fourier_coeffs(com_phase,n_fourier_cutoff);
load('input_insitu_fluctuations.mat')

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
    interference_suite = class_interference_pattern([phase_1(i,:); phase_2(i,:)], [insitu_density_1; insitu_density_2],...
        t_tof);
    condensate_length = interference_suite.condensate_length_Lz;
    grid_dens = interference_suite.output_grid_z; %grid in longitudinal direction
    cloud_widths = interference_suite.compute_density_sigma_t(0, t_tof);  %width of the cloud for imaging input

    rho_tof_full = interference_suite.tof_full_expansion();
    img_rho_tof_full = absorption_imaging(rho_tof_full', grid_dens, cloud_widths);

    z_grid = linspace(-condensate_length/2, condensate_length/2, size(img_rho_tof_full,1));
    x_grid = z_grid;

    density_ripple_no_imaging(i,:) = trapz(interference_suite.output_grid_x, rho_tof_full, 2);
    amp_full = trapz(x_grid, img_rho_tof_full, 2);
    density_ripple(i,:) = amp_full;
    
    count = count+1
end

%Estimating mean density profile - since mean density is estimated, this
%script must run with significantly many number of shots
mean_density = mean(density_ripple, 1); 

%Ignoring the edges
[val, idx_1] = min(abs(z_grid - bulk_start)); %finding start index of the bulk
[val, idx_2] = min(abs(z_grid - bulk_end)); % finding end index of the bulk
cut_coarse_grid_dens = z_grid(idx_1:idx_2); %truncated z grid

cut_density_ripple = density_ripple(:,idx_1:idx_2);
cut_mean_density = mean_density(idx_1:idx_2);
cut_condensate_length = abs(bulk_end-bulk_start);

%Extract common phase
for i = 1:num_samples
    common_suite = class_common_phase_extraction(cut_density_ripple(i,:), cut_mean_density, ...
        cut_coarse_grid_dens, t_tof);
    [cosineCoeffs_u, sineCoeffs_u] = common_suite.extract_com_spectrum(ext_cutoff);
    com_phase_u = common_suite.reconstruct_com_phase(cosineCoeffs_u, sineCoeffs_u);
    output_com_phase(i,:) = com_phase_u;
    output_cosinCoeffs(i,:) = cosineCoeffs_u;
    output_sineCoeffs(i,:) = sineCoeffs_u;
end


for i = 1:10
    figure
    plot(grid_dens, com_phase(i,:))
    hold on
    plot(cut_coarse_grid_dens,output_com_phase(i,:))
end

%save('main_Data/scan_11ms_50nk_with_imaging', 't_tof', 'temp_com', 'temp_rel', 'mean_density','cut_mean_density', 'n_fourier_cutoff', 'ext_cutoff', 'grid_dens','cut_coarse_grid_dens', 'output_com_phase', ...
%   'com_phase', 'rel_phase', 'com_density_fluct', 'rel_density_fluct', 'flag_density_fluct', 'num_samples', 'pixnumz', 'density_ripple', 'cut_density_ripple', 'output_cosinCoeffs', 'output_sineCoeffs',...
%   'in_sineCoeffs', 'in_cosineCoeffs')