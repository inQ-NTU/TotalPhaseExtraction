clear all; close all;
addpath('../../classes/')
addpath('../../plotting_func/')
addpath('results_shaking/scan7436/')
load('scan7436_exp_data_for_phase_extraction_all.mat')

%load density ripple data
density_ripple_all = data_structure_basic.density_profiles_full_stat;
mean_density = data_structure_basic.density_profiles(:,1);

%load parameter
evol_time = data_structure_basic.time;
z = data_structure_basic.z_axis;
t_tof = 11.2e-3;

%set extraction parameter
n_max_fourier = 20;
idx_begin = 20; %The start of untruncated index
idx_end = 64; %The end of untruncated index
z_cut = z(idx_begin:idx_end);
shift = (z_cut(end) - z_cut(1))/2-z_cut(end);
z_cut = z_cut+shift;
dz = abs(z_cut(2)-z_cut(1));
condensate_length = data_structure_basic.L_fit;
mean_density_cut = mean_density(idx_begin:idx_end);

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
        com_ext_suite = class_common_phase_extraction(density_ripple_cut(:,i), mean_density_cut, z_cut, condensate_length, t_tof);
        [output_cosineCoeffs, output_sineCoeffs] = com_ext_suite.extract_com_spectrum_uniform(n_max_fourier);
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
std_cosineCoeffs_all = zeros(length(evol_time), n_max_fourier);
avg_sineCoeffs_all = zeros(length(evol_time), n_max_fourier);
std_sineCoeffs_all = zeros(length(evol_time), n_max_fourier);
%Computing average signal
for i = 1:length(evol_time)
    avg_ext_com_phase_all(i,:) = mean(ext_com_phase_all{i}, 1);
    avg_cosineCoeffs_all(i,:) = mean(ext_cosineCoeffs_all{i}, 1);
    std_cosineCoeffs_all(i,:) = std(ext_cosineCoeffs_all{i},1);
    avg_sineCoeffs_all(i,:) = mean(ext_sineCoeffs_all{i}, 1);
    std_sineCoeffs_all(i,:) = std(ext_sineCoeffs_all{i}, 1);
end

%analyze thermal state
thermal_corr = class_1d_correlation(ext_com_phase_all{1});

%cosine correlation
g1_corr = thermal_corr.g1_corr();
g1_midcut = g1_corr(floor(length(z_cut)/2),:);


%fit the cosine correlation
z_cut = z_cut*1e6;
fitfun = fittype( @(a,b,c,x)a*exp(-abs(x-c)/b));
[fitted_curve,gof] = fit(z_cut(5:end-5)',g1_midcut(5:end-5)',fitfun);
coeffvals = coeffvalues(fitted_curve);

%extract the temperature
hbar = 1.05457e-34; %hbar
m = 86.909*1.66054e-27; %mass of Rb-87 in kg
kb = 1.380649e-23; %Boltzmann constant
n0 = mean(mean_density_cut); %mean density around 70 per microns
%n0 = 70e6;
lambdaT = coeffvals(2)*1e-6; % thermal coherence length in meter

extracted_temperature = 2*(hbar^2)*n0/(m*kb*lambdaT);
confidence_interval = confint(fitted_curve);
lambdaT_uncertainty = (confidence_interval(2,2) - confidence_interval(1,2))*1e-6/2;
temperature_uncertainty = (extracted_temperature/lambdaT)*lambdaT_uncertainty;

%plotting
idx_1 = 9;
idx_2 = 15;
idx_3 = 21;
refined_z = linspace(z_cut(1), z_cut(end), 100);
evol_time = evol_time*1e3;

h = tight_subplot(3,2, [0.1, 0.12], [0.1, 0.05], [0.1, 0.05]);

axes(h(1))
plot(z_cut, ext_com_phase_all{1}, 'Color', [0.7, 0.7,0.7])
hold on
plot(z_cut, avg_ext_com_phase_all(1,:), 'Color','red', 'LineWidth',1.1)
ylim([-4,4])
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
ylim([0,1.1])
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
xlim([-22,22])
ylim([-5,5])
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

axes(h(4))
plot(z_cut, avg_ext_com_phase_all(idx_1,:),'Color', 'red', 'LineWidth',1.1)
hold on
plot(z_cut, avg_ext_com_phase_all(idx_2,:),'-.',  'Color', 'red', 'LineWidth',1.1)
plot(z_cut, avg_ext_com_phase_all(idx_3,:),'--', 'Color', 'red', 'LineWidth',1.1)
xlim([-22,22])
ylim([-3.2,3.2])
xlabel('$z\; (\rm \mu m)$','Interpreter','latex')
ylabel('$\langle\phi_+(z)\rangle$', 'Interpreter','latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);


axes(h(5))
refined_evol_time = linspace(evol_time(1), evol_time(end), 100);
plot(evol_time, avg_cosineCoeffs_all(:,1), 'o', 'Color', 'red', 'LineWidth',1.1)
hold on
plot(refined_evol_time, spline(evol_time, avg_cosineCoeffs_all(:,1), refined_evol_time), 'Color','red', 'LineWidth',1.1)
plot(evol_time, avg_cosineCoeffs_all(:,2)','.-', 'Color', 'blue', 'LineWidth',1.1, 'MarkerSize', 8)
yline(0)
ylim([-3,3])
xline(0, 'LineStyle', '--')
xlim([-35,120])
ylabel('$\langle\rm{Re}[A_k]\rangle$', 'Interpreter','latex')
xlabel('$t\; (\rm ms)$', 'Interpreter','latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{e}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

axes(h(6))
plot(evol_time, avg_sineCoeffs_all(:,1), 'o', 'Color','red', 'LineWidth',1.1)
hold on
plot(refined_evol_time, spline(evol_time, avg_sineCoeffs_all(:,1), refined_evol_time), 'Color','red', 'LineWidth',1.1)
plot(evol_time, avg_sineCoeffs_all(:,2)', '.-', 'Color', 'blue', 'LineWidth',1.1, 'MarkerSize',8)
%plot(evol_time, mode_sine_all(:,3), 'x-', 'Color', 'black')
yline(0)
ylim([-3,3])
xline(0, 'LineStyle', '--')
legend('$k = 2\pi/L$', '', '$k = 4\pi/L$','', '', '', 'Interpreter', 'latex', 'FontSize', 11)
xlim([-35,120])
ylabel('$\langle\rm{Im}[A_k]\rangle$', 'Interpreter','latex')
xlabel('$t\; (\rm ms)$', 'Interpreter','latex')
box on
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{f}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

set(h, 'FontName', 'Times', 'FontSize', 16)