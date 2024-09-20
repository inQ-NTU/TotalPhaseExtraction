clear all; close all;
addpath('../../classes')

load('density_ripple_data_scan5100.mat')
t_tof = 15.6e-3;
ext_cutoff = 20;

for i = 1:length(evol_time)
    mean_density_t = mean(density_ripple_all{i}, 1);
    mean_density{i} = mean_density_t;
    density_ripple = density_ripple_all{i};
    for j = 1:num_of_repetitions
        %com_suite = class_common_phase_extraction(density_ripple(j,:), mean_density_n0, z_axis, t_tof);
        com_suite = class_common_phase_extraction(density_ripple(j,:), mean_density{i}, z_axis, t_tof);
        [cosineCoeffs, sineCoeffs] = com_suite.extract_com_spectrum(ext_cutoff);
        com_phase(j,:) = com_suite.reconstruct_com_phase(cosineCoeffs, sineCoeffs);
    end
    com_phase_all{i} = com_phase;
    mean_com_phase{i} = mean(com_phase, 1);
end
z_axis = com_suite.z_grid;
%plot(z_axis*1e6, mean_density{1}*1e-6, 'o-', 'LineWidth',1.05)
%hold on
%plot(z_axis*1e6, mean_density{2}*1e-6, 'x-', 'LineWidth',1.05)
%plot(z_axis*1e6, mean_density{4}*1e-6, '.-', 'LineWidth',1.05)
%plot(z_axis*1e6, mean_density{5}*1e-6, '^-', 'LineWidth',1.05)
%plot(z_axis*1e6, mean_density{9}*1e-6, 's-', 'LineWidth',1.05)
%legend('$t = -1.6\; \rm ms$', '$t = 0\; \rm ms$', '$t = 5\; \rm ms$', '$t = 20\; \rm ms$', '$t = 30\; \rm ms$', 'Interpreter', 'latex',...
%    'FontSize', 12)
%ylabel('$n_0 \; \rm (\mu m^{-1})$', 'Interpreter','latex', 'FontSize',18)
%xlabel('$z\; \rm (\mu m)$', 'Interpreter','latex', 'FontSize',18)
%save('ext_com_phase_data_scan5100', 'com_phase_all', 'mean_com_phase', 'mean_density', 'density_ripple_all', 't_tof', 'ext_cutoff', ...
%    'z_axis', 'scan_no', 'evol_time','num_of_repetitions')