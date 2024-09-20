clear all; close all;
addpath('../plotting_func/')
load('sample_thermometry_50nk_100microns.mat')
g1_L1 = g1;
fit_curve_L1 = fitted_curve;
corr_theory_L1 = corr_theory;
T1 = extracted_temperature;
z_grid_L1 = z_grid;

load('sample_thermometry_50nk_300microns.mat')
g1_L2 = g1;
fit_curve_L2 = fitted_curve;
corr_theory_L2 = corr_theory;
T2 = extracted_temperature;
z_grid_L2 = z_grid;

load('sample_thermometry_50nk_500microns.mat')
g1_L3 = g1;
fit_curve_L3 = fitted_curve;
corr_theory_L3 = corr_theory;
T3 = extracted_temperature;
z_grid_L3 = z_grid;

f = tight_subplot(1,3,[.06 .06],[.2 .1],[.1 .1]);

axes(f(1))
plot(z_grid_L1, g1_L1, 'o', 'MarkerSize', 4)
hold on
plot(z_grid_L1, fit_curve_L1(z_grid_L1), 'LineWidth',1.05);
plot(z_grid_L1, corr_theory_L1, 'LineWidth',1.05, 'LineStyle','--', 'Color', 'black')
%legend('Sampled (N = 10^4)', "Fitted (T = "+num2str(T1*1e9, 3)+" nK)",...
%    'Theory (Discrete)', 'FontSize', 10)
%title('$T = 20 \; \rm nK$', 'Interpreter','latex')
%legend box off
ylabel('$C_+(z)$', 'Interpreter', 'latex')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylim([0,1])
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);

axes(f(2))
plot(z_grid_L2, g1_L2, 'o', 'MarkerSize',4)
hold on
plot(z_grid_L2, fit_curve_L2(z_grid_L2), 'LineWidth',1.05);
plot(z_grid_L2, corr_theory_L2, 'LineWidth',1.05, 'LineStyle','--', 'Color', 'black')
%legend('Sampled', "Fitted (T = "+num2str(T2*1e9, 3)+ " nK)",...
%    'Theory (Discrete)', 'FontSize', 10)
%title('$T = 40 \; \rm nK$', 'Interpreter','latex')
%legend box off
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylim([0,1])
yticks([])
xlim([-150,150])
xticks([-150,0,150])
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);

axes(f(3))
plot(z_grid_L3, g1_L3, 'o',  'MarkerSize', 4)
hold on
plot(z_grid_L3, fit_curve_L3(z_grid_L3), 'LineWidth',1.05);
plot(z_grid_L3, corr_theory_L3, 'LineWidth',1.05, 'LineStyle','--', 'Color','black')
%legend('Sampled (N = 10^4)', "Fitted (T = "+num2str(T3*1e9, 3)+ " nK)",...
%   'Theory (Discrete)', 'FontSize', 10)
%title('$T = 60\; \rm nK$', 'Interpreter','latex')
%legend box off
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylim([0,1])
yticks([])
xticks([-250,0,250])
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);

set(f, 'FontName', 'Times', 'FontSize', 14)