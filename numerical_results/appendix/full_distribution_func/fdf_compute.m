clear all;
close all;
addpath('../../main_Data/')

load('scan_11ms_16ms_50nk_80microns_10ksamples.mat')

integration_length_l = 60e-6;

pixel_size = abs(cut_grid_dens(2)-cut_grid_dens(1));
midpoint = floor(length(cut_grid_dens)/2);

idx_begin = midpoint - floor(integration_length_l/(2*pixel_size));
idx_end = midpoint + floor(integration_length_l/(2*pixel_size));

com_phase = com_phase(:,idx_1:idx_2);

for i = 1:num_samples
    truncated_in_com_phase = com_phase(i,idx_begin:idx_end);
    truncated_out_com_phase_t1 = out_com_phase_t1(i,idx_begin:idx_end);
    truncated_out_com_phase_t2 = out_com_phase_t2(i,idx_begin:idx_end);

    exp_phase_in = arrayfun(@(x) exp(1j*x), truncated_in_com_phase);
    exp_phase_out_t1 = arrayfun(@(x) exp(1j*x), truncated_out_com_phase_t1);
    exp_phase_out_t2 = arrayfun(@(x) exp(1j*x), truncated_out_com_phase_t2);
    
    integrated_contrast_in(i) = abs(trapz(pixel_size, exp_phase_in)).^2;
    integrated_contrast_out_t1(i) = abs(trapz(pixel_size, exp_phase_out_t1)).^2;
    integrated_contrast_out_t2(i) = abs(trapz(pixel_size, exp_phase_out_t2)).^2;

end

integrated_contrast_in = integrated_contrast_in/mean(integrated_contrast_in);
integrated_contrast_out_t1 = integrated_contrast_out_t1/mean(integrated_contrast_out_t1);
integrated_contrast_out_t2 = integrated_contrast_out_t2/mean(integrated_contrast_out_t2);

figure
h1 = histogram(integrated_contrast_in, 'BinWidth',0.1, 'Normalization','pdf', 'DisplayStyle','stairs', 'LineWidth',1.1);
hold on
histogram(integrated_contrast_out_t1, 'BinWidth',0.1, 'Normalization','pdf', 'DisplayStyle', 'stairs', 'LineWidth',1.1)
histogram(integrated_contrast_out_t2, 'BinWidth',0.1, 'Normalization','pdf', 'DisplayStyle','stairs', 'LineWidth',1.1)

save('fdf_60microns.mat', 'integrated_contrast_in', 'integrated_contrast_out_t1', 'integrated_contrast_out_t2')