clear all;close all;
addpath('Data_new/')
addpath('../plotting_func/')
addpath('../classes')

load('scan_20ms_50nk_no_imaging.mat')
input_com_phase_20ms = com_phase;
output_com_phase_20ms = output_com_phase;

load('scan_11ms_50nk_no_imaging.mat')
input_com_phase_11ms = com_phase;
output_com_phase_11ms = output_com_phase;

%load('scan_11ms_50nk_no_imaging_no_DF.mat')
%output_com_phase_11ms_no_DF = output_com_phase;


corr_in = class_1d_correlation(input_com_phase_11ms);
corr_out_11ms = class_1d_correlation(output_com_phase_11ms);
%corr_out_11ms_no_DF = class_1d_correlation(output_com_phase_11ms_no_DF);
corr_out_20ms = class_1d_correlation(output_com_phase_20ms);

g1_in = corr_in.g1_corr();
g1_out_11ms = corr_out_11ms.g1_corr();
g1_out_20ms = corr_out_20ms.g1_corr();

mode_energy_in = diag(corr_in.fourier_correlation());
mode_energy_out_11ms = diag(corr_out_11ms.fourier_correlation());
%mode_energy_out_11ms_no_DF = diag(corr_out_11ms_no_DF.fourier_correlation());
mode_energy_out_20ms = diag(corr_out_20ms. fourier_correlation());

out_pixnum = length(cut_grid_dens);
in_pixnum = length(grid_dens);

f = tight_subplot(1,2, [0.15, 0.12], [0.2, 0.1], [0.12, 0.05]);

axes(f(1))
plot(grid_dens*1e6, g1_in(floor(in_pixnum/2),:),'Color','Black', 'LineStyle','--', 'LineWidth',1.1)
hold on
plot(cut_grid_dens*1e6, g1_out_11ms(floor(out_pixnum/2),:), 'Color', 'red', 'LineWidth',1.1)
plot(cut_grid_dens*1e6,g1_out_20ms(floor(out_pixnum/2),:),'Color', 'Blue', 'LineWidth',1.1)
ylim([0,1])
ax = gca;
ax.LineWidth = 1.1;
ylabel('$C_+(z)$', 'Interpreter', 'latex')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);

axes(f(2))
plot(1:5,mode_energy_in(2:6), 'o--', 'Color', 'Black', 'LineWidth',1.1)
hold on
plot(1:5,mode_energy_out_20ms(2:6), '^-b', 'LineWidth',1.1)
plot(1:5,mode_energy_out_11ms(2:6), 'x-r', 'LineWidth',1.1, 'MarkerSize', 8)
%plot(1:5,mode_energy_out_11ms_no_DF(2:6), 's:r', 'LineWidth',1.5)
%legend('Input', '$t_{\rm tof} = 20\; \rm ms$','$t_{\rm tof} = 11\; \rm ms$', '$t_{\rm tof} = 11\; \rm{ms}, \; \delta n_+ = 0$',...
%    'FontSize',10, 'Interpreter', 'latex')
%legend box off
xlim([1,5])
ax = gca;
ax.LineWidth = 1.1;
ylabel('$\langle |A_p|^2\rangle$', 'Interpreter', 'latex')
xlabel('Mode index $p$', 'Interpreter', 'latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.85]);
ylim([0,0.8])
g= axes('Position',[.72 .55 .2 .3]);
box on
plot(3:8,mode_energy_in(4:9), 'o-.','Color', 'Black', 'LineWidth',1.1)
hold on
plot(3:8,mode_energy_out_20ms(4:9), '^-b', 'LineWidth',1.1)
plot(3:8,mode_energy_out_11ms(4:9), 'x-r', 'MarkerSize', 8, 'LineWidth',1.1)
%plot(3:8,mode_energy_out_11ms_no_DF(4:9), 's:r', 'LineWidth',1.5)
ax = gca;
ax.LineWidth = 1.1;
xticks([3,4,5,6,7,8])
yticks([0,0.02,0.04,0.06])
xlim([3,8])
set(g, 'FontName', 'Times', 'FontSize',14)

set(f, 'FontName', 'Times', 'FontSize', 18)
