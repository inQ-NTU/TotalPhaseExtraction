clear all;
close all;
addpath('../../../plotting_func')

load('fdf_20microns.mat')
fdf_in_20microns = integrated_contrast_in;
fdf_out_20microns_t1 = integrated_contrast_out_t1;
fdf_out_20microns_t2 = integrated_contrast_out_t2;

load('fdf_60microns.mat')
fdf_in_60microns = integrated_contrast_in;
fdf_out_60microns_t1 = integrated_contrast_out_t1;
fdf_out_60microns_t2 = integrated_contrast_out_t2;

f = tight_subplot(1,2, [0.08, 0.1], [0.25, 0.08], [0.1, 0.08]);

axes(f(1))
histogram(fdf_in_20microns, 'Normalization','pdf' ,'BinWidth',0.1, 'Normalization', 'pdf', 'DisplayStyle','stairs', 'EdgeColor','black',...
    'LineWidth', 1.05, 'LineStyle', '-.')
hold on
histogram(fdf_out_20microns_t1, 'Normalization','pdf', 'BinWidth',0.1, 'Normalization', 'pdf', 'DisplayStyle','stairs', 'EdgeColor', 'red', ...
    'LineWidth', 1.05)
histogram(fdf_out_20microns_t2, 'Normalization','pdf', 'BinWidth',0.1, 'Normalization', 'pdf', 'DisplayStyle','stairs', 'EdgeColor', 'blue', ...
    'LineWidth', 1.05)
xlabel('$\xi_+$', 'Interpreter','latex')
ylabel('$P(\xi_+)$', 'Interpreter', 'latex')
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.9,0.82]);
ax = gca;
ax.LineWidth = 1.1;
ylim([0,2.5])
xlim([0,3])
%ylim([0,1.25])
%ylim([0,2.1])

axes(f(2))
histogram(fdf_in_60microns, 'Normalization','pdf' ,'BinWidth',0.1, 'Normalization', 'pdf', 'DisplayStyle','stairs', 'EdgeColor','black', 'LineWidth', 1.05,...
    'LineStyle','-.')
hold on
histogram(fdf_out_60microns_t1, 'Normalization','pdf', 'BinWidth',0.1, 'Normalization', 'pdf', 'DisplayStyle','stairs', 'EdgeColor', 'red', ...
    'LineWidth', 1.05)
histogram(fdf_out_60microns_t2, 'Normalization','pdf', 'BinWidth',0.1, 'Normalization', 'pdf', 'DisplayStyle','stairs', 'EdgeColor', 'blue', ...
    'LineWidth', 1.05)
xlabel('$\xi_+$', 'Interpreter','latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.9,0.82]);
ax = gca;
ax.LineWidth = 1.1;
%ylim([0,1.25])
%yticks([])
%ylim([0,2.1])
xlim([0,3])
ylim([0,0.8])
yticks([0,0.4,0.8])
%ylabel('$P(\xi_+)$', 'Interpreter', 'latex')

set(f, 'FontName', 'Times', 'FontSize', 16)
