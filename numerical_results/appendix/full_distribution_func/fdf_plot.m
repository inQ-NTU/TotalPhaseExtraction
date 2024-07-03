clear all;
close all;
addpath('../../../plotting_func')

load('fdf_20microns.mat')
fdf_in_20microns = integrated_contrast_in;
fdf_out_20microns = integrated_contrast_out;

load('fdf_30microns.mat')
fdf_in_30microns = integrated_contrast_in;
fdf_out_30microns = integrated_contrast_out;

load('fdf_60microns.mat')
fdf_in_60microns = integrated_contrast_in;
fdf_out_60microns = integrated_contrast_out;

f = tight_subplot(1,3, [0.08, 0.05], [0.25, 0.08], [0.1, 0.08]);

axes(f(1))
histogram(fdf_in_20microns, 'Normalization','pdf' ,'BinWidth',0.2)
hold on
histogram(fdf_out_20microns, 'Normalization','pdf', 'BinWidth',0.2)
xlabel('$\xi_+$', 'Interpreter','latex')
ylabel('$P(\xi_+)$', 'Interpreter', 'latex')
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.9,0.82]);
ax = gca;
ax.LineWidth = 1.1;
xlim([0,3])
text(1.5,0.9,'$l = 20\; \rm \mu m$', 'Interpreter', 'latex', 'FontSize',14)
ylim([0,1.25])
axes(f(2))
histogram(fdf_in_30microns, 'Normalization','pdf', 'BinWidth',0.2)
hold on
histogram(fdf_out_30microns, 'Normalization','pdf', 'BinWidth',0.2)
xlabel('$\xi_+$', 'Interpreter','latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.9,0.82]);
ax = gca;
ax.LineWidth = 1.1;
ylim([0,1.25])
xlim([0,3])
yticks([])
text(1.6,0.9,'$l = 30\; \rm \mu m$', 'Interpreter', 'latex', 'FontSize',14)

axes(f(3))
histogram(fdf_in_60microns, 'Normalization','pdf', 'BinWidth',0.2)
hold on
histogram(fdf_out_60microns, 'Normalization','pdf', 'BinWidth',0.2)
xlabel('$\xi_+$', 'Interpreter','latex')
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.9,0.82]);
ax = gca;
ax.LineWidth = 1.1;
ylim([0,1.25])
yticks([])
xlim([0,3])
text(1.6,0.9,'$l = 60\; \rm \mu m$', 'Interpreter', 'latex', 'FontSize',14)

set(f, 'FontName', 'Times', 'FontSize', 16)
