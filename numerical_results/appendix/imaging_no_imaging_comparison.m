clear all; close all;
addpath('../../plotting_func/')
addpath('../../classes')
addpath('../main_Data/')

load('scan_11ms_16ms_50nK_80microns.mat')
out_com_phase_no_imaging = out_com_phase_t1;
density_ripple_no_imaging = density_ripple_t1*1e-6;
cosineCoeffs_no_imaging = out_cosineCoeffs_t1;
sineCoeffs_no_imaging = out_sineCoeffs_t1;
spectrum_no_imaging = cosineCoeffs_no_imaging.^2+sineCoeffs_no_imaging.^2;
mean_spectrum_no_imaging = mean(spectrum_no_imaging, 1);
input_spectrum = in_cosineCoeffs.^2 +in_sineCoeffs.^2;
mean_input_spectrum = mean(input_spectrum, 1);

corr_input = class_1d_correlation(com_phase);
g1_input = corr_input.g1_corr();
g1_input = g1_input(floor(length(grid_dens)/2), :);

corr_no_imaging = class_1d_correlation(out_com_phase_no_imaging);
g1_no_imaging = corr_no_imaging.g1_corr();
g1_no_imaging = g1_no_imaging(floor(length(cut_grid_dens)/2),:);

load('scan_11ms_50nk_with_imaging.mat')
spectrum = output_cosinCoeffs.^2 + output_sineCoeffs.^2;
mean_spectrum = mean(spectrum, 1);
corr = class_1d_correlation(output_com_phase);
g1 = corr.g1_corr();
g1 = g1(floor(length(cut_coarse_grid_dens)/2),:);

%Plotting
grid_dens = grid_dens*1e6;
%coarse_grid_dens = coarse_grid_dens*1e6;
cut_coarse_grid_dens = cut_coarse_grid_dens*1e6;
cut_grid_dens = cut_grid_dens*1e6;
chosen_idx = 3;

f = tight_subplot(2,2, [0.15, 0.12], [0.2, 0.1], [0.12, 0.05]);

axes(f(1))
plot(grid_dens, density_ripple_no_imaging(chosen_idx,:), 'Color','blue', 'LineWidth',1.05)
yline(150, 'LineStyle','-.', 'LineWidth',1.05);
ylabel('$n_{\rm tof}\; (\rm \mu m^{-1})$', 'Interpreter','latex')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylim([0,300])
yticks([0,100,200,300])
hold on 
yyaxis right
plot(cut_coarse_grid_dens, cut_density_ripple(chosen_idx,:), 'Color', 'red', 'LineWidth',1.05)
xline(-40, 'LineStyle','--', 'LineWidth',1.05)
xline(40, 'LineStyle','--', 'LineWidth',1.05)
ax = gca;
ax.YAxis(2).Exponent = -2;
ax.YAxis(2).Color = 'red';
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.15,0.8]);


axes(f(2))
plot(cut_grid_dens, out_com_phase_no_imaging(chosen_idx,:), 'Color', 'blue', 'LineWidth',1.05)
hold on
plot(cut_coarse_grid_dens, output_com_phase(chosen_idx,:), 'Color', 'red', 'LineWidth',1.05)
plot(grid_dens, com_phase(chosen_idx,:), 'Color','black', 'LineStyle','-.','LineWidth',1.05)
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$\phi_+(z)$', 'Interpreter','latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);


axes(f(3))
plot(cut_grid_dens, g1_no_imaging, 'Color','blue', 'LineWidth',1.05)
hold on
plot(cut_coarse_grid_dens, g1, 'Color', 'red', 'LineWidth',1.05)
plot(grid_dens, g1_input, 'Color','black', 'LineStyle','-.', 'LineWidth',1.05)
xline(-40, 'LineStyle','--', 'LineWidth',1.05)
xline(40, 'LineStyle','--', 'LineWidth',1.05)
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$C_+(z)$', 'Interpreter','latex')
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.15,0.8]);


axes(f(4))
plot(1:5, mean_spectrum(1:5), '^-r', 'LineWidth',1.05)
hold on
plot(1:5, mean_spectrum_no_imaging(1:5), 'x-b', 'LineWidth', 1.05)
plot(1:5, mean_input_spectrum(1:5), '-.', 'Color','black', 'LineWidth',1.05)
legend('With Imaging', 'Without Imaging', 'Input', 'FontSize', 12)
xlabel('Mode index $p$', 'Interpreter', 'latex')
ylabel('$\langle |A_p|^2\rangle$', 'Interpreter','latex')
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);

set(f, 'FontName', 'Times', 'FontSize', 14)