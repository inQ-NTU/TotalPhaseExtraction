close all; clear all;
addpath('../../classes/')
addpath('../../plotting_func/')

load('ext_com_phase_data_scan5100.mat')

%extract the temperature
z_axis = z_axis*1e6;
hbar = 1.05457e-34; %hbar
m = 86.909*1.66054e-27; %mass of Rb-87 in kg
kb = 1.380649e-23; %Boltzmann constant
fitfun = fittype( @(a,b,c,x)a*exp(-abs(x-c)/b));
fit_idx = 7;

for i = 1:length(evol_time)

    %Computing correlation function
    corr_suite = class_1d_correlation(com_phase_all{i});
    g1 = corr_suite.g1_corr();
    cov = corr_suite.covariance_matrix();
    cov_matrix{i} = cov;
    g1 = g1(ceil(length(z_axis)/2),:);
    g1_all{i} = g1;

    %Fit thermal coherence length
    [fitted_curve, gof] = fit(z_axis(fit_idx:end-fit_idx)', g1(fit_idx:end-fit_idx)', fitfun, 'StartPoint', [1,10,0]);
    coeffvals = coeffvalues(fitted_curve);
    confidence_interval = confint(fitted_curve);
    coh_length = coeffvals(2);
    lambdaT(i) = coh_length*1e-6;
    DlambdaT_upper(i) = (confidence_interval(2,2)-coh_length)*1e-6;
    DlambdaT_lower(i)  = (coh_length-confidence_interval(1,2))*1e-6;
    
    %Extract temperature
    n0 = mean(mean_density{i});
    temperature_T(i) = (hbar^2)*n0/(2*m*kb*lambdaT(i));
    Dtemperature_T_upper(i) = (temperature_T(i)/lambdaT(i))*DlambdaT_upper(i);
    Dtemperature_T_lower(i) = (temperature_T(i)/lambdaT(i))*DlambdaT_lower(i);
end

%Plotting start
h = tight_subplot(2,2, [0.12, 0.12], [0.15, 0.1], [0.1, 0.05]);

idx_0 = 1;
idx_1 = 4;
idx_2 = 6;
idx_3 = 8;

axes(h(1))
plot(z_axis, com_phase_all{idx_0}, 'Color', [.7 .7 .7])
hold on
plot(z_axis, mean_com_phase{idx_0}, 'Color', 'red', 'LineWidth', 1.08)
xlabel('$z\; (\rm \mu m)$', 'Interpreter', 'latex')
ylabel('$\phi_+(z)$', 'Interpreter', 'latex')
ylim([-4,4])
xlim([-36,36])

axes(h(2))
plot(z_axis, mean_com_phase{idx_0}, 'Color', 'red', 'LineWidth', 1.08)
hold on
plot(z_axis, mean_com_phase{idx_1}, 'Color', 'blue', 'LineWidth', 1.08, 'LineStyle','--')
plot(z_axis, mean_com_phase{idx_2}, 'Color', "#77AC30", 'LineWidth', 1.08, 'LineStyle','-.')
plot(z_axis, mean_com_phase{idx_3}, 'Color', "black", 'LineWidth', 1.5, 'LineStyle',':')
xlabel('$z\; (\rm \mu m)$', 'Interpreter', 'latex')
ylabel('$\langle \phi_+(z)\rangle$', 'Interpreter', 'latex')
%ylim([-8e-16, 8e-16])

axes(h(3))
plot(z_axis, g1_all{idx_0}, 'Color', 'red', 'LineWidth', 1.08)
hold on
plot(z_axis, g1_all{idx_1}, 'Color','blue', 'LineWidth', 1.08, 'LineStyle','--')
plot(z_axis, g1_all{idx_2}, 'Color',"#77AC30", 'LineWidth', 1.08, 'LineStyle','-.')
plot(z_axis, g1_all{idx_3}, 'Color', "black", 'LineWidth', 1.5, 'LineStyle',':')
xlabel('$z\; (\rm \mu m)$', 'Interpreter', 'latex')
ylabel('$C_+(z)$', 'Interpreter','latex')
xline(z_axis(fit_idx), 'LineStyle','--')
xline(z_axis(end-fit_idx), 'LineStyle','--')
legend('$t = 0\; \rm ms$', '$t = 6\; \rm ms$', '$t = 12\; \rm ms$', '$t = 18\; \rm ms$', 'Interpreter', 'latex', 'FontSize', 10)
legend box off

axes(h(4))
errorbar(evol_time, temperature_T*1e9, Dtemperature_T_lower*1e9, Dtemperature_T_upper*1e9, 'ob')
%xlim([-3,66])
%xlim([-3,34])
%ylim([14,35])
%ylim([0,20])
ylim([15,30])
xlim([-3,19])
xlabel('$t\; (\rm ms)$', 'Interpreter','latex')
ylabel('$T_+\; (\rm nK)$', 'Interpreter','latex')
xline(0, 'LineStyle','--')

set(h, 'FontName', 'Times', 'FontSize', 14)