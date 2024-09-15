clear all; 
close all;

addpath('../../classes')

load('ext_com_phase_bootstrap_data_scan5100.mat')

time_index = 0;
for i = 1:4
    bootstrap_datasets = bootstrap_com_phase_all{i};
    time_index = time_index+1
    bootstrap_index = 0;
    tic
    for j = 1:length(bootstrap_datasets)
        bootstrap_index = bootstrap_index + 1
        ext_com_phases = bootstrap_datasets{j}';
        corr_suite = class_1d_correlation(ext_com_phases);
        g4_disconnected{j} = corr_suite.wick_four_point_correlation();
        g4_full{j} = corr_suite.correlation_func(4);
        toc
        tic
    end
    g4_disconnected_all{i} = g4_disconnected;
    g4_full_all{i} = g4_full;
end

save('g4_data_bootstrap_scan5100_evol_time_idx_1_to_4.mat', 'bootstrap_com_phase_all', 'com_phase_all', 'g4_full_all', 'g4_disconnected_all', 'evol_time', 'z_axis')

