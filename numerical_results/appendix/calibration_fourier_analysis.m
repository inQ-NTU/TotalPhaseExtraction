clear all;close all;
addpath('../../classes')
addpath('../../plotting_func')
addpath('../Data')
load('scan_11ms_calibrate_LE.mat')
com_phase_in = com_phase;
com_phase_out_LE = output_com_phase;
density_ripple_example_1 = density_ripple;
output_com_phase_example_1 = output_com_phase;

load('scan_11ms_calibrate_LE_RP.mat');
com_phase_out_LE_and_RP = output_com_phase;
density_ripple_example_2 = density_ripple;
output_com_phase_example_2 = output_com_phase;

load('scan_11ms_calibrate_LE_RP_DF.mat')
com_phase_out_LE_and_DF_and_RP = output_com_phase;

corr_in = class_1d_correlation(com_phase_in);
corr_LE = class_1d_correlation(com_phase_out_LE);
corr_LE_and_RP = class_1d_correlation(com_phase_out_LE_and_RP);
corr_LE_and_DF_and_RP = class_1d_correlation(com_phase_out_LE_and_DF_and_RP);

fourier_in = diag(corr_in.fourier_correlation());
fourier_LE = diag(corr_LE.fourier_correlation());
fourier_LE_and_RP = diag(corr_LE_and_RP.fourier_correlation());
fourier_LE_and_DF_and_RP = diag(corr_LE_and_DF_and_RP.fourier_correlation());

idx = 999;
grid_dens = grid_dens*1e6;
cut_grid_dens = cut_grid_dens*1e6;
f = tight_subplot(1,3, [0.08, 0.1], [0.3, 0.1], [0.12, 0.12]);
axes(f(1))
plot(grid_dens, density_ripple_example_1(idx,:)*1e-6,'Color','blue', 'LineWidth',1.1)
hold on
plot(grid_dens, density_ripple_example_2(idx,:)*1e-6,'Color','red', 'LineWidth',1.1)
yline(2*mean_density, 'LineWidth',1.1, 'LineStyle','-.','Color','black')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$n_{\rm tof}(\rm \mu m^{-1})$', 'Interpreter','latex')
ylim([50,250])
yticks([50,150,250])
yline(150, '-.', 'LineWidth',1.1)
xline(-40, '--', 'LineWidth',1.1)
xline(40, '--', 'LineWidth',1.1)
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);

axes(f(2))
plot(grid_dens, com_phase(idx,:), 'Color','Black', 'LineStyle', '-.', 'LineWidth',1.1)
hold on
plot(cut_grid_dens, output_com_phase_example_1(idx,:), 'Color','blue', 'LineWidth',1.1)
plot(cut_grid_dens, output_com_phase_example_2(idx,:), 'Color','red', 'LineWidth',1.1)
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$\phi_+(z)$', 'Interpreter','latex')
xline(-40, '--', 'LineWidth',1.1)
xline(40, '--', 'LineWidth',1.1)
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);

ax = gca;
ax.LineWidth = 1.1;

axes(f(3))
%plot(fourier_in(2:5), 'o-', 'Color', 'Black', 'LineStyle','--', 'LineWidth',1.1)
plot(1:4, fourier_LE(2:5), 'o-','Color','Blue','LineWidth',1.1)
hold on
plot(1:4, fourier_LE_and_RP(2:5),'x-', 'Color','Black', 'LineWidth',1.1,'MarkerSize', 8)
plot(1:4, fourier_LE_and_DF_and_RP(2:5),'^-', 'Color','red', 'LineWidth',1.1)
xlim([1,4])
legend('LE', 'LE + RP', 'LE + RP + DF', 'FontSize', 10)
xlabel('Mode index $p$', 'Interpreter','latex')
ylabel('$\langle |A_p|^2\rangle$', 'Interpreter','latex')
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);
ax = gca;
ax.LineWidth = 1.1;
ax.YAxis.Exponent = -2;
ylim([0,0.09])
xticks([1,2,3,4])
yticks([0,0.03,0.06,0.09])

set(f, 'FontName', 'Times', 'FontSize', 16)
