close all;clear all;
addpath('../classes')
addpath('../plotting_func/')
addpath('data/')

%load 20 ms TOF data
load('scan_20ms_5nk.mat')
fid_20ms_5nk = output_fidelity;

load('scan_20ms_10nk.mat')
fid_20ms_10nk = output_fidelity;

load('scan_20ms_15nk.mat')
fid_20ms_15nk = output_fidelity;

load('scan_20ms_20nk.mat')
fid_20ms_20nk = output_fidelity;
out_com_phase = output_com_phase;
input_com_phase = com_phase;

load('scan_20ms_25nk.mat')
fid_20ms_25nk = output_fidelity;

load('scan_20ms_30nk.mat')
fid_20ms_30nk = output_fidelity;


%load different expansion time data
%5 nK dataset
load('scan_11ms_5nK.mat')
fid_11ms_5nk = output_fidelity;

load('scan_15ms_5nK.mat')
fid_15ms_5nk = output_fidelity;

load('scan_30ms_5nK.mat')
fid_30ms_5nk = output_fidelity;

load('scan_40ms_5nK.mat')
fid_40ms_5nk = output_fidelity;

%20 nK dataset
load('scan_11ms_20nK.mat')
fid_11ms_20nk = output_fidelity;

load('scan_15ms_20nK.mat')
fid_15ms_20nk = output_fidelity;

load('scan_30ms_20nK.mat')
fid_30ms_20nk = output_fidelity;

load('scan_40ms_20nK.mat')
fid_40ms_20nk = output_fidelity;

%median fidelity (interquartile range) vs expansion time 
t_tof = [11,15,20,30,40];
%median fidelity
median_fid_5nk = [median(fid_11ms_5nk), median(fid_15ms_5nk), ...
    median(fid_20ms_5nk), median(fid_30ms_5nk), median(fid_40ms_5nk)];
median_fid_20nk = [median(fid_11ms_20nk), median(fid_15ms_20nk), ...
    median(fid_20ms_20nk), median(fid_30ms_20nk), median(fid_40ms_20nk)];

%interquartile range
upper_error_5nK = [quantile(fid_11ms_5nk,0.75) - median(fid_11ms_5nk), ...
    quantile(fid_15ms_5nk,0.75) - median(fid_15ms_5nk), ...
    quantile(fid_20ms_5nk,0.75) - median(fid_20ms_5nk), ...
    quantile(fid_30ms_5nk,0.75) - median(fid_30ms_5nk), ...
    quantile(fid_40ms_5nk,0.75) - median(fid_40ms_5nk), ...
    ];

upper_error_20nK = [quantile(fid_11ms_20nk,0.75) - median(fid_11ms_20nk), ...
    quantile(fid_15ms_20nk,0.75) - median(fid_15ms_20nk), ...
    quantile(fid_20ms_20nk,0.75) - median(fid_20ms_20nk), ...
    quantile(fid_30ms_20nk,0.75) - median(fid_30ms_20nk), ...
    quantile(fid_40ms_20nk,0.75) - median(fid_40ms_20nk), ...
    ];

lower_error_5nK = [-quantile(fid_11ms_5nk,0.25) + median(fid_11ms_5nk), ...
    -quantile(fid_15ms_5nk,0.25) + median(fid_15ms_5nk), ...
    -quantile(fid_20ms_5nk,0.25) + median(fid_20ms_5nk), ...
    -quantile(fid_30ms_5nk,0.25) + median(fid_30ms_5nk), ...
    -quantile(fid_40ms_5nk,0.25) + median(fid_40ms_5nk), ...
    ];

lower_error_20nk = [-quantile(fid_11ms_20nk,0.25) + median(fid_11ms_20nk), ...
    -quantile(fid_15ms_20nk,0.25) + median(fid_15ms_20nk), ...
    -quantile(fid_20ms_20nk,0.25) + median(fid_20ms_20nk), ...
    -quantile(fid_30ms_20nk,0.25) + median(fid_30ms_20nk), ...
    -quantile(fid_40ms_20nk,0.25) + median(fid_40ms_20nk), ...
    ];

temperature_list = [5,10,15,20,25,30];
median_fid_20ms = [median(fid_20ms_5nk), median(fid_20ms_10nk), median(fid_20ms_15nk), ...
    median(fid_20ms_20nk), median(fid_20ms_25nk), median(fid_20ms_30nk)];
upper_error_20ms = [quantile(fid_20ms_5nk, 0.75) - median(fid_20ms_5nk),...
    quantile(fid_20ms_10nk,0.75) - median(fid_20ms_10nk), ...
    quantile(fid_20ms_15nk,0.75) - median(fid_20ms_15nk), ...
    quantile(fid_20ms_20nk,0.75) - median(fid_20ms_20nk), ...
    quantile(fid_20ms_25nk,0.75) - median(fid_20ms_25nk), ...
    quantile(fid_20ms_30nk,0.75) - median(fid_20ms_30nk), ...
    ];

lower_error_20ms = [-quantile(fid_20ms_5nk, 0.25) + median(fid_20ms_5nk),...
    -quantile(fid_20ms_10nk,0.25) + median(fid_20ms_10nk), ...
    -quantile(fid_20ms_15nk,0.25) + median(fid_20ms_15nk), ...
    -quantile(fid_20ms_20nk,0.25) + median(fid_20ms_20nk), ...
    -quantile(fid_20ms_25nk,0.25) + median(fid_20ms_25nk), ...
    -quantile(fid_20ms_30nk,0.25) + median(fid_20ms_30nk), ...
    ];

%analysis correlation
corr_out = class_1d_correlation(out_com_phase);
corr_in = class_1d_correlation(input_com_phase);

g1_out = corr_out.g1_corr();
g1_in = corr_in.g1_corr();
out_pixnum = length(cut_z_grid);
in_pixnum = length(z_grid);

f = tight_subplot(2,2, [0.15, 0.12], [0.15, 0.1], [0.12, 0.05]);

axes(f(1))
histogram(fid_20ms_5nk, 'Normalization','probability')
xlim([0.85,1])
ylabel('$P(F)$', 'Interpreter','latex')
xlabel('$F$', 'Interpreter','latex')
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);
box on
ax = gca;
ax.LineWidth = 1.1;


axes(f(2))
errorbar(t_tof, median_fid_5nk,lower_error_5nK,upper_error_5nK, 'o--', 'Color', 'blue','LineWidth',1.1)
hold on
errorbar(t_tof, median_fid_20nk, lower_error_20nk, upper_error_20nK,  '^--', 'Color','blue', 'LineWidth',1.1)
ylim([0.6,1.05])
xlim([9,42])
yline(1, '-.')
xlabel('$t_{\rm TOF}\; (\rm ms)$', 'Interpreter', 'latex')
ylabel('$\overline{F}$', 'Interpreter','latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);
box on
ax = gca;
ax.LineWidth = 1.1;

axes(f(3))
errorbar(temperature_list, median_fid_20ms,lower_error_20ms, upper_error_20ms, '--square', 'Color','blue', 'LineWidth',1.1)
xlabel('$T_+ \; (\rm nK)$', 'Interpreter','latex')
ylabel('$\overline{F}$', 'Interpreter','latex')
xlim([4,31])
ylim([0.7,1.05])
xticks([5,10,15,20, 25,30])
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);
box on
ax = gca;
ax.LineWidth = 1.1;

axes(f(4))
plot(z_grid*1e6, g1_in(in_pixnum/2,:), 'Color', 'black', 'LineWidth',1.1, 'LineStyle','-.')
hold on
plot(cut_z_grid*1e6, g1_out(out_pixnum/2,:), 'Color','black', 'LineWidth',1.1)
xlabel('$z\; \rm (\mu m)$', 'Interpreter','latex')
ylabel('$C(z)$', 'Interpreter','latex')
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.1,0.8]);
box on
ax = gca;
ax.LineWidth = 1.1;

set(f, 'FontName', 'Times', 'FontSize', 16)

