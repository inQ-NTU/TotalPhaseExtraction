clear all; close all;
addpath('Data_new/')
addpath('../plotting_func/')
load('scan_11ms_50nk_no_imaging.mat')

chosen_idx = 456;
grid_dens = grid_dens*1e6;
cut_grid_dens = cut_grid_dens*1e6;

f = tight_subplot(2,1,[.08 .08],[.15 .1],[.2 .1]);
axes(f(1))
plot(grid_dens, density_ripple(chosen_idx,:)*1e-6, 'LineWidth',1.2, 'Color','red')
yline(2*mean_density*1e-6, '-.', 'LineWidth',1.2, 'LineStyle','-.')
xline(-40, '--', 'LineWidth',1.2)
xline(40, '--', 'LineWidth',1.2)
xticks([])
yticks([0,150,300])
ylim([0,350])
xlim([-50,50])
ylabel('$n_{\rm tof}\;(\rm \mu m^{-1})$', 'Interpreter','latex')
ax = gca;
ax.LineWidth = 1.2;

axes(f(2))
plot(cut_grid_dens,output_com_phase(chosen_idx,:), 'Color','Red', 'LineWidth',1.2)
hold on
plot(grid_dens, com_phase(chosen_idx,:), 'Color','Black', 'LineStyle','-.',...
    'LineWidth',1.2)
xline(-40, '--', 'LineWidth',1.2)
xline(40, '--', 'LineWidth',1.2)
ylabel('$\phi_+(z)$', 'Interpreter','latex')
xlabel('$z\;(\rm \mu m)$', 'Interpreter','latex')
ax = gca;
ax.LineWidth = 1.2;

set(f, 'FontName', 'Times', 'FontSize', 20)
