clear all;
close all;
addpath('../../Data/')

load('scan_11ms_50nk_no_imaging.mat')

integration_length_l = 30e-6;

pixel_size = abs(grid_dens(2)-grid_dens(1));
integration_pixnum = floor(integration_length_l/pixel_size);
midpoint_out = floor(length(cut_grid_dens)/2);
midpoint_in = floor(length(grid_dens)/2);

idx_begin_out = midpoint_out - floor(integration_pixnum/2);
idx_end_out = midpoint_out + floor(integration_pixnum/2);

idx_begin_in = midpoint_in - floor(integration_pixnum/2);
idx_end_in = midpoint_in + floor(integration_pixnum/2);


for i = 1:num_samples
    cut_in_com_phase = com_phase(i,idx_begin_in:idx_end_in);
    cut_out_com_phase = com_phase(i, idx_begin_out:idx_end_out);

    exp_phase_in = arrayfun(@(x) exp(1j*x), cut_in_com_phase);
    exp_phase_out = arrayfun(@(x) exp(1j*x), cut_out_com_phase);
    
    integrated_contrast_in(i) = abs(trapz(pixel_size, exp_phase_in)).^2;
    integrated_contrast_out(i) = abs(trapz(pixel_size, exp_phase_out)).^2;

end

integrated_contrast_in = integrated_contrast_in/mean(integrated_contrast_in);
integrated_contrast_out = integrated_contrast_out/mean(integrated_contrast_out);

histogram(integrated_contrast_in)
hold on
histogram(integrated_contrast_out)

save('fdf_30microns.mat', 'integrated_contrast_in', 'integrated_contrast_out')
