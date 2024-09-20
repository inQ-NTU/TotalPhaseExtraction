clear all; close all;
addpath('../../classes/')
addpath('../../plotting_func/')
addpath('results_shaking/scan7436/')
load('scan7436_exp_data_for_phase_extraction_all.mat')

thermal_state_idx = 2;
%load density ripple data
density_ripple_all = data_structure_basic.density_profiles_full_stat;
mean_density = data_structure_basic.density_profiles(:,thermal_state_idx);

%load parameter
evol_time = data_structure_basic.time;
z = data_structure_basic.z_axis;
t_tof = 11.2e-3;

%set extraction parameter
n_max_fourier = 40;
idx_begin = 25; %The start of untruncated index
idx_end = 60;
z_cut = z(idx_begin:idx_end);
%shift = (z_cut(end) - z_cut(1))/2-z_cut(end);
shift = z_cut(1);
z_cut = z_cut-shift;
%z_cut = z_cut+shift;
dz = abs(z_cut(2)-z_cut(1));

mean_density_cut = mean_density(idx_begin:idx_end);
n0 = mean(mean_density_cut); %mean density around 64 per microns

ext_com_phase_all = {};
ext_cosineCoeffs_all = {};
ext_sineCoeffs_all = {};

for j = 1:length(evol_time)
    density_ripple = density_ripple_all{j};
    num_samples = size(density_ripple, 2);
    %Truncate the boundary
    density_ripple_cut = density_ripple(idx_begin:idx_end,:);
    ext_com_phase = zeros(num_samples, length(z_cut));
    ext_cosineCoeffs = zeros(num_samples, n_max_fourier);
    ext_sineCoeffs = zeros(num_samples, n_max_fourier);
    for i = 1:num_samples
        com_ext_suite = class_common_phase_extraction(density_ripple_cut(:,i), mean_density_cut, z_cut, t_tof,...
            'Neumann');
        [output_cosineCoeffs, output_sineCoeffs] = com_ext_suite.extract_com_spectrum(n_max_fourier);
        ext_cosineCoeffs(i,:) = output_cosineCoeffs;
        ext_sineCoeffs(i,:) = output_sineCoeffs;
        ext_com_phase(i,:) = com_ext_suite.reconstruct_com_phase(output_cosineCoeffs, output_sineCoeffs, z_cut);
    end
    ext_com_phase_all{j} = ext_com_phase;
    ext_cosineCoeffs_all{j} = ext_cosineCoeffs;
    ext_sineCoeffs_all{j} = ext_sineCoeffs;
end

avg_ext_com_phase_all = zeros(length(evol_time), length(z_cut));
avg_cosineCoeffs_all = zeros(length(evol_time), n_max_fourier);
avg_sineCoeffs_all = zeros(length(evol_time), n_max_fourier);

%Computing average signal
for i = 1:length(evol_time)
    avg_ext_com_phase_all(i,:) = mean(ext_com_phase_all{i}, 1);
    avg_cosineCoeffs_all(i,:) = mean(ext_cosineCoeffs_all{i}, 1);
    avg_sineCoeffs_all(i,:) = mean(ext_sineCoeffs_all{i}, 1);
end

%analyze thermal state
thermal_corr = class_1d_correlation(ext_com_phase_all{thermal_state_idx});

%cosine correlation
g1_corr = thermal_corr.g1_corr();
g1_midcut = g1_corr(floor(length(z_cut)/2),:);

%fit the cosine correlation
z_cut = z_cut*1e6;
fitfun = fittype( @(b,c,x)exp(-abs(x-c)/b));
[fitted_curve,gof] = fit(z_cut(10:end-10)',g1_midcut(10:end-10)',fitfun, 'StartPoint', [10, 20]);
coeffvals = coeffvalues(fitted_curve);

%extract the temperature
hbar = 1.05457e-34; %hbar
m = 86.909*1.66054e-27; %mass of Rb-87 in kg
kb = 1.380649e-23; %Boltzmann constant
n0 = mean(mean_density_cut); %mean density around 64 per microns
lambdaT = coeffvals(1)*1e-6; % thermal coherence length in meter

extracted_temperature = 2*(hbar^2)*n0/(m*kb*lambdaT);
confidence_interval = confint(fitted_curve);
lambdaT_uncertainty = (confidence_interval(2,2) - confidence_interval(1,2))*1e-6/2;
temperature_uncertainty = (extracted_temperature/lambdaT)*lambdaT_uncertainty;

%plotting
idx_1 = 9;
idx_2 = 15;
idx_3 = 21;
evol_time = evol_time*1e3;
refined_z = linspace(z_cut(1), z_cut(end), 100);

figure
h = tight_subplot(3,2, [0.12, 0.12], [0.15, 0.1], [0.1, 0.05]);

axes(h(1))
plot(z_cut, ext_com_phase_all{thermal_state_idx}, 'Color', [0.7, 0.7,0.7])
hold on
plot(z_cut, avg_ext_com_phase_all(thermal_state_idx,:), 'Color','red', 'LineWidth',1.1)
%ylim([-4,4])
xlim([0,48])
%xticks([])
xlabel('$z\; (\rm \mu m)$', 'Interpreter', 'latex')
ylabel('$\phi_+(z)$', 'Interpreter','latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

axes(h(2))
plot(z_cut, g1_midcut,'x', 'Color','black', 'LineWidth',1.1)
hold on
plot(refined_z, fitted_curve(refined_z), 'Color','black', 'LineWidth',1.1)
yticks([0,0.5,1])
xlim([0,48])
xline(z_cut(15))
xline(z_cut(end-15))
%xticks([])
%ylim([0,1.1])
ylabel('$C_+(z)$', 'Interpreter','latex')
xlabel('$z\; (\rm \mu m)$', 'Interpreter', 'latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

axes(h(3))
plot(z_cut, ext_com_phase_all{idx_1}, 'Color', [0.7,0.7,0.7])
hold on
plot(z_cut, avg_ext_com_phase_all(idx_1,:),'Color','red', 'LineWidth',1.1)
xlabel('$z\; (\rm \mu m)$','Interpreter','latex')
ylabel('$\phi_+(z)$', 'Interpreter','latex')
xlim([0,48])
%ylim([-5,5])
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

axes(h(4))
plot(z_cut, avg_ext_com_phase_all(idx_1,:),'Color', 'red', 'LineWidth',1.1)
hold on
plot(z_cut, avg_ext_com_phase_all(idx_2,:),'-.',  'Color', 'red', 'LineWidth',1.1)
plot(z_cut, avg_ext_com_phase_all(idx_3,:),'--', 'Color', 'red', 'LineWidth',1.1)
xlim([0,48])
%ylim([-3.2,3.2])
xlabel('$z\; (\rm \mu m)$','Interpreter','latex')
ylabel('$\langle\phi_+(z)\rangle$', 'Interpreter','latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

axes(h(5))
refined_evol_time = linspace(evol_time(1), evol_time(end), 100);
plot(evol_time, avg_cosineCoeffs_all(:,2), 'o', 'Color', 'red', 'LineWidth',1.1, 'MarkerSize', 4)
hold on
plot(refined_evol_time, spline(evol_time, avg_cosineCoeffs_all(:,2), refined_evol_time), 'Color','red', 'LineWidth',1.1)
plot(evol_time, avg_cosineCoeffs_all(:,1)','.-', 'Color', 'blue', 'LineWidth',1.1, 'MarkerSize', 8)
yline(0)
%ylim([-3,3])
yticks([-3,3])
xline(0, 'LineStyle', '--')
xline(evol_time(idx_1), 'LineStyle', '--')
xline(evol_time(idx_2), 'LineStyle','--')
xline(evol_time(idx_3), 'LineStyle','--')
xlim([-35,120])
ylabel('$\langle\rm{Re}[A_k]\rangle$', 'Interpreter','latex')
xlabel('$t\; (\rm ms)$', 'Interpreter','latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{e}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

axes(h(6))
plot(evol_time, avg_sineCoeffs_all(:,1), 'o', 'Color','red', 'LineWidth',1.1, 'MarkerSize', 4)
hold on
plot(refined_evol_time, spline(evol_time, avg_sineCoeffs_all(:,1), refined_evol_time), 'Color','red', 'LineWidth',1.1)
plot(evol_time, avg_sineCoeffs_all(:,2)', '.-', 'Color', 'blue', 'LineWidth',1.1, 'MarkerSize',8)
%plot(evol_time, mode_sine_all(:,3), 'x-', 'Color', 'black')
yline(0)
ylim([-3,3])
yticks([-3,3])
xline(0, 'LineStyle', '--')
xline(evol_time(idx_1), 'LineStyle', '--')
xline(evol_time(idx_2), 'LineStyle','--')
xline(evol_time(idx_3), 'LineStyle','--')
legend('$k = 2\pi/L$', '', '$k = 4\pi/L$','', '', '', 'Interpreter', 'latex', 'FontSize', 11)
legend('boxoff')
xlim([-35,120])
ylabel('$\langle\rm{Im}[A_k]\rangle$', 'Interpreter','latex')
xlabel('$t\; (\rm ms)$', 'Interpreter','latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{f}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

set(h, 'FontName', 'Times', 'FontSize', 16)