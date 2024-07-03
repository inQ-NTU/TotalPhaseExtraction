clear all; close all;
addpath('../../plotting_func/')
addpath('../../classes')
addpath('../Data/')

load('input_comparison_imaging.mat')
input_com_phase = com_phase;

load('scan_11ms_50nk_no_imaging.mat')
output_com_phase_no_imaging = output_com_phase;
output_density_ripple_no_imaging = density_ripple*1e-6;
mean_density_ori = mean_density*1e-6; 

load('scan_11ms_50nk_with_imaging.mat')
output_com_phase_with_imaging = output_com_phase;
output_density_ripple_with_imaging = density_ripple;

chosen_idx = 786;
condensate_length = 100e-6;
coarse_grid_dens = linspace(-condensate_length/2, condensate_length/2, size(output_density_ripple_with_imaging,2));

%Correlations
corr_in = class_1d_correlation(input_com_phase);
corr_out_no_imaging = class_1d_correlation(output_com_phase_no_imaging);
corr_out_with_imaging = class_1d_correlation(output_com_phase_with_imaging);

g1_in = corr_in.g1_corr();
g1_out_no_imaging = corr_out_no_imaging.g1_corr();
g1_out_with_imaging = corr_out_with_imaging.g1_corr();

fourier_in = diag(corr_in.fourier_correlation());
fourier_out_no_imaging = diag(corr_out_no_imaging.fourier_correlation());
fourier_out_with_imaging = diag(corr_out_with_imaging.fourier_correlation());

in_pixnum = length(grid_dens);
out_pixnum_no_imaging = length(cut_grid_dens);
out_pixnum_with_imaging = length(cut_coarse_grid_dens);

%Plotting
grid_dens = grid_dens*1e6;
coarse_grid_dens = coarse_grid_dens*1e6;
cut_coarse_grid_dens = cut_coarse_grid_dens*1e6;
cut_grid_dens = cut_grid_dens*1e6;

f = tight_subplot(2,2, [0.15, 0.12], [0.2, 0.1], [0.12, 0.05]);

axes(f(1))
plot(grid_dens, output_density_ripple_no_imaging(chosen_idx,:), 'Color','red', 'LineWidth',1.1)
yline(2*mean_density_ori, 'Color', 'Black', 'LineStyle','-.', 'LineWidth',1.1)
xline(-40, '--', 'LineWidth',1.1)
xline(40, '--', 'LineWidth',1.1)
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$n_{\rm tof}\; (\rm \mu m^{-1})$', 'Interpreter','latex')
ylim([0,350])
hold on
yyaxis right
plot(coarse_grid_dens, output_density_ripple_with_imaging(chosen_idx,:), 'Color', 'blue', 'LineWidth',1.1)
ax = gca;
ax.LineWidth = 1.1;
ax.YAxis(2).Color = 'blue';
ax.YAxis(2).Exponent = -2;
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);



axes(f(2))
plot(grid_dens, input_com_phase(chosen_idx,:), 'Color','Black', 'LineStyle','-.', 'LineWidth',1.1)
hold on
plot(cut_grid_dens, output_com_phase_no_imaging(chosen_idx,:), 'Color', 'red', 'LineWidth',1.1)
plot(cut_coarse_grid_dens, output_com_phase_with_imaging(chosen_idx,:), 'Color','blue', 'LineWidth',1.1)
xline(-40, '--', 'LineWidth',1.1)
xline(40, '--', 'LineWidth',1.1)
ax = gca;
ax.LineWidth = 1.1;
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$\phi_+(z)$', 'Interpreter','latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);



axes(f(3))
plot(grid_dens, g1_in(:,floor(in_pixnum/2)), 'Color','black', 'LineStyle','-.', 'LineWidth',1.1);
hold on
plot(cut_grid_dens, g1_out_no_imaging(:,floor(out_pixnum_no_imaging/2)), 'Color','red', 'LineWidth',1.1)
plot(cut_coarse_grid_dens, g1_out_with_imaging(:,floor(out_pixnum_with_imaging/2)), 'Color','blue', 'LineWidth',1.1)
ax = gca; 
ax.LineWidth = 1.1;
ylim([0,1])
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylabel('$C_+(z)$', 'Interpreter','latex')
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);


axes(f(4))
plot(1:5, fourier_in(2:6), 'o-','Color', 'black', 'LineWidth',1.1, 'LineStyle','--')
hold on
plot(1:5, fourier_out_no_imaging(2:6), 'x-', 'Color','red', 'MarkerSize',8, 'LineWidth',1.1)
plot(1:5, fourier_out_with_imaging(2:6), '^-', 'Color','blue', 'LineWidth',1.1)
legend('Input', 'Without Imaging', 'With Imaging','FontSize',12)
xticks([1,2,3,4,5])
xlim([1,5])
xlabel('Mode index $p$', 'Interpreter', 'latex')
ylabel('$\langle |A_p|^2\rangle$', 'Interpreter','latex')
ax = gca;
ax.LineWidth = 1.1;
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);
ylim([0,0.82])
yticks([0,0.2,0.4,0.6,0.8])
set(f, 'FontName', 'Times', 'FontSize', 16)