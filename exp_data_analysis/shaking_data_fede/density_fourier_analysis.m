clear all; close all;
addpath('../../classes/')
addpath('../../plotting_func/')
addpath('results_shaking/scan7436/')
load('scan7436_exp_data_for_phase_extraction_all.mat')

%load time and z grid
evol_time = data_structure_basic.time;
z = data_structure_basic.z_axis;

thermal_state_idx = 1;
n_max_fourier = 20;
idx_begin = 20; %The start of untruncated index
idx_end = 64; %The end of untruncated index
z_cut = z(idx_begin:idx_end);
shift = z_cut(1);
z_cut = z_cut-shift;

%load density ripple data
density_ripple_all = data_structure_basic.density_profiles_full_stat;
mean_density_all = data_structure_basic.density_profiles(idx_begin:idx_end,:);
mean_density_thermal = data_structure_basic.density_profiles(idx_begin:idx_end,thermal_state_idx);

for i = 1:length(evol_time)
    density_perturbation(i,:) = (mean_density_all(:,i) - mean_density_thermal)./mean_density_thermal;
end

%Cosine mode analysis of density perturbation
cosineCoeffs_cosineSeries = zeros(length(evol_time), n_max_fourier);
cosineCoeffs_fourierSeries = zeros(length(evol_time), n_max_fourier);
sineCoeffs_fourierSeries = zeros(length(evol_time), n_max_fourier);
for i = 1:length(evol_time)
    cosineCoeffs_cosineSeries(i,:) = cosine_decomposition(density_perturbation(i,:), z_cut, n_max_fourier);
    [cosineCoeffs_fourierSeries(i,:), sineCoeffs_fourierSeries(i,:)] = ...
        fourier_decomposition(density_perturbation(i,:), z_cut, n_max_fourier);
end

evol_time = evol_time*1e3;
idx_1 = 9;
idx_2 = 15;
idx_3 = 21;

h = tight_subplot(2,2, [0.1, 0.1], [0.15, 0.1], [0.1, 0.05]);

axes(h(1))
plot(evol_time, cosineCoeffs_cosineSeries(:,2), 'o-', 'MarkerSize', 4, 'Color', 'red', 'LineWidth',1.1)
hold on
plot(evol_time, cosineCoeffs_cosineSeries(:,4), '.-', 'Color', 'blue', 'LineWidth',1.1)
ylim([-0.2,0.2])
xline(0, 'LineStyle', '--')
xline(evol_time(idx_1), 'LineStyle', '--')
xline(evol_time(idx_2), 'LineStyle','--')
xline(evol_time(idx_3), 'LineStyle','--')
yline(0)
xticks([])
ylabel('$A_k$', 'Interpreter', 'latex')
legend('$k = 2\pi/L$', '$k = 4\pi/L$', 'Interpreter', 'latex', 'FontSize', 12)

axes(h(2))
plot(evol_time, cosineCoeffs_cosineSeries(:,1), 'o-', 'MarkerSize', 4, 'Color', 'red', 'LineWidth',1.1)
hold on
plot(evol_time, cosineCoeffs_cosineSeries(:,3), '.-', 'Color', 'blue', 'LineWidth',1.1)
ylim([-0.2, 0.2])
xline(0, 'LineStyle', '--')
xline(evol_time(idx_1), 'LineStyle', '--')
xline(evol_time(idx_2), 'LineStyle','--')
xline(evol_time(idx_3), 'LineStyle','--')
yline(0)
xticks([])
ylabel('$A_k$', 'Interpreter', 'latex')
legend('$k = \pi/L$', '$k = 3\pi/L$', 'Interpreter', 'latex', 'FontSize', 12)

axes(h(3))
plot(evol_time, cosineCoeffs_fourierSeries(:,1), 'o-', 'MarkerSize', 4, 'Color', 'red', 'LineWidth',1.1)
hold on
plot(evol_time, cosineCoeffs_fourierSeries(:,2), '.-', 'Color', 'blue', 'LineWidth',1.1)
ylim([-0.2,0.2])
xline(0, 'LineStyle', '--')
xline(evol_time(idx_1), 'LineStyle', '--')
xline(evol_time(idx_2), 'LineStyle','--')
xline(evol_time(idx_3), 'LineStyle','--')
yline(0)
ylabel('$\rm{Re}(A_k)$', 'Interpreter', 'latex')
legend('$k = 2\pi/L$', '$k = 4\pi/L$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('$t\; (\rm ms)$', 'Interpreter', 'latex')

axes(h(4))
plot(evol_time, sineCoeffs_fourierSeries(:,1), 'o-', 'MarkerSize', 4, 'Color', 'red', 'LineWidth',1.1)
hold on
plot(evol_time, sineCoeffs_fourierSeries(:,2), '.-', 'Color', 'blue', 'LineWidth',1.1)
ylim([-0.2,0.2])
xline(0, 'LineStyle', '--')
xline(evol_time(idx_1), 'LineStyle', '--')
xline(evol_time(idx_2), 'LineStyle','--')
xline(evol_time(idx_3), 'LineStyle','--')
yline(0)
ylabel('$\rm{Im}(A_k)$', 'Interpreter', 'latex')
legend('$k = 2\pi/L$', '$k = 4\pi/L$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('$t\; (\rm ms)$', 'Interpreter','latex')

set(h, 'FontName', 'Times', 'FontSize', 16)