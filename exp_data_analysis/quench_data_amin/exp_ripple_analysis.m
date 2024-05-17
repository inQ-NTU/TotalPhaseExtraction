clear all
close all
addpath('../../classes/')

load('uncoupled_to_coupled_DW_box.mat')
ext_cutoff = 20;
condensate_length = (z_um(end)-z_um(1))*1e-6;
t_tof = 15.6e-3;

%-2 ms data (beginning of the cooling), uncoupled double well (J = )
density_ripple = density_arr{:,2}; 
mean_density = mean_density_per_um{:,2};
phase_data = phase_arr{2,:};


%Truncating boundary and converting grid to SI
pix_cut = 6; %Cut 5 pixels on each boundary
density_ripple = density_ripple(:,pix_cut:end-pix_cut);
mean_density = mean_density(pix_cut:end-pix_cut);

z_grid = z_um.*1e-6;
z_grid_cut = z_um(pix_cut:end-pix_cut).*1e-6;

%Computing scaled density ripple
num_samples = size(density_ripple, 1);
pix_num = size(density_ripple, 2);
scaled_density_ripple = zeros(num_samples, pix_num);
for i = 1:num_samples
    scaled_density_ripple(i,:) = 1-(density_ripple(i,:)./mean_density);
end


%Extracting common phase
ext_com_phase_all_shots = zeros(num_samples, length(z_grid));
ext_cosineCoeffs_all_shots = zeros(num_samples, ext_cutoff);
ext_sineCoeffs_all_shots = zeros(num_samples, ext_cutoff);

for i = 1:num_samples
    com_phase_suite = class_common_phase_spectrum(scaled_density_ripple(i,:), z_grid_cut, t_tof, condensate_length);
    [ext_cosineCoeffs, ext_sineCoeffs] = com_phase_suite.extract_com_spectrum(ext_cutoff);
    ext_com_phase = com_phase_suite.extract_com_profile(z_grid);
    ext_com_phase_all_shots(i,:) = ext_com_phase;
    ext_cosineCoeffs_all_shots(i,:) = ext_cosineCoeffs;
    ext_sineCoeffs_all_shots(i,:) = ext_sineCoeffs;
end

%save('ext_com_phase_data_uncoupled_batch2.mat', 'ext_com_phase_all_shots', 'ext_cosineCoeffs_all_shots', 'ext_sineCoeffs_all_shots', 'ext_cutoff', 'z_grid')







