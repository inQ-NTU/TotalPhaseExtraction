clear all;close all;
%addpath('main_Data/')
addpath('main_Data/')
addpath('../plotting_func/')
addpath('../classes')
load('scan_11ms_16ms_50nk_80microns.mat')
input_com_phase = com_phase(:,idx_1:idx_2);

out_pixnum = length(cut_grid_dens);
corr_in = class_1d_correlation(input_com_phase);
corr_out_t1 = class_1d_correlation(out_com_phase_t1);
corr_out_t2 = class_1d_correlation(out_com_phase_t2);

g1_in = corr_in.g1_corr();
g1_out_t1 = corr_out_t1.g1_corr();
g1_out_t2 = corr_out_t2.g1_corr();

%Extracting only the midcut
g1_in_midcut = g1_in(floor(out_pixnum/2),:);
g1_out_t1_midcut = g1_out_t1(floor(out_pixnum/2),:);
g1_out_t2_midcut = g1_out_t2(floor(out_pixnum/2),:);

%Computing mode energy
mode_energy_in = in_cosineCoeffs.^2+in_sineCoeffs.^2;
mode_energy_t1 = out_cosineCoeffs_t1.^2+out_sineCoeffs_t1.^2;                           
mode_energy_t2 = out_cosineCoeffs_t2.^2+out_sineCoeffs_t2.^2;

mean_mode_energy_in = mean(mode_energy_in, 1);
mean_mode_energy_t1 = mean(mode_energy_t1, 1);
mean_mode_energy_t2 = mean(mode_energy_t2, 1);

%Temperature extraction
cut_grid_dens = cut_grid_dens*1e6;
grid_dens = grid_dens*1e6;
fitfun = fittype(@(a,b,c,x)a*exp(-abs(x-c)/b));

%Fitting
fit_start = -20;
fit_end = 20;
[val, fit_idx_1] = min(abs(cut_grid_dens - fit_start)); %finding start index of the bulk
[val, fit_idx_2] = min(abs(cut_grid_dens - fit_end)); % finding end index of the bulk
[fitted_curve_in,gof_in] = fit(cut_grid_dens(fit_idx_1:fit_idx_2)',g1_in_midcut(fit_idx_1:fit_idx_2)',fitfun);
[fitted_curve_t1,gof_t1] = fit(cut_grid_dens(fit_idx_1:fit_idx_2)',g1_out_t1_midcut(fit_idx_1:fit_idx_2)',fitfun);
[fitted_curve_t2,gof_t2] = fit(cut_grid_dens(fit_idx_1:fit_idx_2)',g1_out_t2_midcut(fit_idx_1:fit_idx_2)',fitfun);

coeffvals_in = coeffvalues(fitted_curve_in);
coeffvals_out_t1 = coeffvalues(fitted_curve_t1);
coeffvals_out_t2 = coeffvalues(fitted_curve_t2);

hbar = 1.05457e-34; %hbar
m = 86.909*1.66054e-27; %mass of Rb-87 in kg
kb = 1.380649e-23; %Boltzmann constant 
lambda_t1 = coeffvals_out_t1(2)*1e-6; % thermal coherence length in meter
lambda_t2 = coeffvals_out_t2(2)*1e-6;
lambda_in = coeffvals_in(2)*1e-6;

extracted_temperature_in = (hbar^2)*mean_density/(m*kb*lambda_in);
extracted_temperature_out_t1 = (hbar^2)*mean_density/(m*kb*lambda_t1);
extracted_temperature_out_t2 = (hbar^2)*mean_density/(m*kb*lambda_t2);

%Coarse graining the correlation data to make the plot less crowded
cg_cut_grid_dens = [];
cg_g1_in_midcut = [];
cg_g1_out_t1_midcut = [];
cg_g1_out_t2_midcut = [];
for i = 1:length(cut_grid_dens)
    if rem(i, 2) 
        cg_cut_grid_dens = [cg_cut_grid_dens, cut_grid_dens(i)];
        cg_g1_in_midcut = [cg_g1_in_midcut, g1_in_midcut(i)];
        cg_g1_out_t1_midcut = [cg_g1_out_t1_midcut, g1_out_t1_midcut(i)];
        cg_g1_out_t2_midcut = [cg_g1_out_t2_midcut, g1_out_t2_midcut(i)];
    end
end

f = tight_subplot(1,2,[0.15, 0.12], [0.2, 0.1], [0.12, 0.05]);

axes(f(1))
plot(cg_cut_grid_dens, cg_g1_in_midcut,'o', 'Color', 'black', 'MarkerSize', 3)
hold on
plot(cut_grid_dens, fitted_curve_in(cut_grid_dens), 'Color','Black', 'LineWidth',1.1)
plot(cg_cut_grid_dens, cg_g1_out_t1_midcut, 'xr', 'MarkerSize',5)
plot(cut_grid_dens, fitted_curve_t1(cut_grid_dens), 'Color', 'red', 'LineWidth',1.1)
plot(cg_cut_grid_dens,cg_g1_out_t2_midcut, '^b', 'MarkerSize', 3.5)
plot(cut_grid_dens, fitted_curve_t2(cut_grid_dens), 'Color', 'blue', 'LineWidth',1.1)
xline(cut_grid_dens(fit_idx_1), 'LineStyle','--', 'LineWidth',1.1)
xline(cut_grid_dens(fit_idx_2), 'LineStyle','--', 'LineWidth',1.1)
ylim([0,1.1])
xlim([-40,40])
ax = gca;
ax.LineWidth = 1.1;
ylabel('$C_+(z)$', 'Interpreter', 'latex')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);

axes(f(2))
plot(1:5,mean_mode_energy_in(1:5), 'o-.', 'Color', 'Black', 'LineWidth',1.1)
hold on
plot(1:5,mean_mode_energy_t1(1:5), 'x-r', 'LineWidth',1.1, 'MarkerSize', 7)
plot(1:5,mean_mode_energy_t2(1:5), '^-b', 'LineWidth',1.1, 'MarkerSize', 6)
xlim([1,5])
ax = gca;
ax.LineWidth = 1.1;
ylabel('$\langle |A_p|^2\rangle$', 'Interpreter', 'latex')
xlabel('Mode index $p$', 'Interpreter', 'latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);
g= axes('Position',[.72 .55 .2 .3]);
box on
plot(3:8, mean_mode_energy_in(3:8), 'o-.','Color', 'Black', 'LineWidth',1.1)
hold on
plot(3:8, mean_mode_energy_t1(3:8), 'x-r', 'LineWidth',1.1, 'MarkerSize', 7)
plot(3:8, mean_mode_energy_t2(3:8), '^-b', 'MarkerSize', 5, 'LineWidth',1.1)
ax = gca;
ax.LineWidth = 1.1;
xticks([3,4,5,6,7,8])
xlim([3,8])
set(g, 'FontName', 'Times', 'FontSize',14)

set(f, 'FontName', 'Times', 'FontSize', 18)