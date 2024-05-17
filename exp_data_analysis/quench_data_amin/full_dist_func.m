clear all
close all
load('ext_com_phase_data_uncoupled_batch1.mat')
condensate_length = (z_grid(end) - z_grid(1))*1e6;
sublength = 50;
num_samples = size(ext_com_phase_all_shots,1);
dim = size(ext_com_phase_all_shots, 2);

pixel_size = condensate_length/(dim-1);

mid_point = floor(dim/2);
end_idx_1 = floor(mid_point-sublength*(dim-1)/(2*condensate_length)); 
end_idx_2 = floor(mid_point+sublength*(dim-1)/(2*condensate_length));

integrated_contrast = zeros(1,num_samples);

for i = 1:num_samples
    cut_phase = ext_com_phase_all_shots(i,end_idx_1:end_idx_2);
    exp_phase = arrayfun(@(x) exp(1j*x), cut_phase);
    integrated_contrast(i) = trapz(pixel_size,exp_phase);
end

integrated_contrast = abs(integrated_contrast).^2;
integrated_contrast = integrated_contrast./mean(integrated_contrast);