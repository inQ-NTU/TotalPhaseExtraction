close all;
clear all;

addpath('../../plotting_func/')
addpath('../../classes')
addpath('../main_Data/')

load('scan_11ms_16ms_50nk_LE.mat')
density_ripple_LE = density_ripple_t1;
com_phase_LE = out_com_phase_t1;
cosine_LE = out_cosineCoeffs_t1;
sine_LE = out_sineCoeffs_t1;
mean_spectrum_LE = mean(cosine_LE.^2+sine_LE.^2, 1);

load('scan_11ms_16ms_50nk_LE_RP.mat')
density_ripple_LE_RP = density_ripple_t1;
com_phase_LE_RP = out_com_phase_t1;
cosine_LE_RP = out_cosineCoeffs_t1;
sine_LE_RP = out_sineCoeffs_t1;
mean_spectrum_LE_RP = mean(cosine_LE_RP.^2+sine_LE_RP.^2, 1);

load('scan_11ms_16ms_50nk_LE_RP_DF.mat')
density_ripple_LE_RP_DF = density_ripple_t1;
com_phase_LE_RP_DF = out_com_phase_t1;
cosine_LE_RP_DF = out_cosineCoeffs_t1;
sine_LE_RP_DF = out_sineCoeffs_t1;
mean_spectrum_LE_RP_DF = mean(cosine_LE_RP_DF.^2+sine_LE_RP_DF.^2, 1);

idx = 121;
grid_dens = grid_dens*1e6;
cut_grid_dens = cut_grid_dens*1e6;
f = tight_subplot(1,3, [0.08, 0.1], [0.3, 0.1], [0.12, 0.12]);
axes(f(1))
plot(grid_dens, density_ripple_LE(idx,:)*1e-6,'Color','black', 'LineWidth',1.05)
hold on
plot(grid_dens, density_ripple_LE_RP(idx,:)*1e-6,'Color','blue', 'LineWidth',1.05)
plot(grid_dens, density_ripple_LE_RP_DF(idx,:)*1e-6,'Color','red', 'LineWidth',1.05)
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$n_{\rm tof}(\rm \mu m^{-1})$', 'Interpreter','latex')
ylim([50,250])
%ylim([50,250])
% yticks([50,150,250])
yline(150, '-.')
xline(-40, '--')
xline(40, '--')
ax = gca;
ax.LineWidth = 1.05;
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);


axes(f(2))
plot(cut_grid_dens, com_phase_LE(idx,:), 'Color','black', 'LineWidth',1.05)
hold on
plot(cut_grid_dens, com_phase_LE_RP(idx,:), 'Color','blue', 'LineWidth',1.05)
plot(cut_grid_dens, com_phase_LE_RP_DF(idx,:),'Color','red', 'LineWidth',1.05)
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$\phi_+(z)$', 'Interpreter','latex')
xline(-40, '--', 'LineWidth',1.05)
xline(40, '--', 'LineWidth',1.05)
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);

ax = gca;
ax.LineWidth = 1.05;

axes(f(3))
plot(1:5, mean_spectrum_LE(1:5), 'o-','Color','Black', 'LineWidth',1.05)
hold on
plot(1:5, mean_spectrum_LE_RP(1:5),'x-', 'Color','Blue','MarkerSize', 8, 'LineWidth',1.05)
plot(1:5, mean_spectrum_LE_RP_DF(1:5),'^-', 'Color','red', 'LineWidth',1.05)
xlim([1,4])
legend('LE', 'LE + RP', 'LE + RP + DF', 'FontSize', 8)
xlabel('Mode index $p$', 'Interpreter','latex')
ylabel('$\langle |A_p|^2\rangle$', 'Interpreter','latex')
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.82]);
ax = gca;
ax.LineWidth = 1.05;
xticks([1,2,3,4,5])

set(f, 'FontName', 'Times', 'FontSize', 14)



