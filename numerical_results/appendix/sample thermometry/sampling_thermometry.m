close all; clear all;
addpath('../classes/')
addpath('../plotting_func/')

%set parameter
T = 50e-9;
mean_density = 75e6;
condensate_length = 500e-6;
pixnumz = condensate_length*1e6;
num_samples = 1000;
max_fourier = floor(pixnumz/2);
z_grid = linspace(-condensate_length/2, condensate_length/2, pixnumz);
scaling_factor = sqrt(2);

%begin sampling
sampling_suite = class_bogoliubov_sampling(T, mean_density, scaling_factor, condensate_length);
[phase_sample, density_fluct_sample] = sampling_suite.generate_fluct_samples(max_fourier, pixnumz, num_samples);

corr_theory = sampling_suite.corr_function(max_fourier, pixnumz);
corr_suite = class_1d_correlation(phase_sample);

g1 = corr_suite.g1_corr();
g1 = g1(pixnumz/2,:);

%fit the cosine correlation
z_grid = z_grid*1e6;
%Fitting
fit_start = -20;
fit_end = 20;
[~, fit_idx_1] = min(abs(z_grid - fit_start)); %finding start index of the bulk
[~, fit_idx_2] = min(abs(z_grid - fit_end)); % finding end index of the bulk

fitfun = fittype(@(a,b,c,x)a*exp(-a*abs(x-c)/b));
[fitted_curve,gof] = fit(z_grid(fit_idx_1:fit_idx_2)',g1(fit_idx_1:fit_idx_2)',fitfun);

coeffvals = coeffvalues(fitted_curve);

%extract the temperature
hbar = 1.05457e-34; %hbar
m = 86.909*1.66054e-27; %mass of Rb-87 in kg
kb = 1.380649e-23; %Boltzmann constant
lambdaT = coeffvals(2)*1e-6; % thermal coherence length in meter

extracted_temperature = (hbar^2)*mean_density/(m*kb*lambdaT);

save('sample_thermometry_50nk_500microns.mat', "extracted_temperature", 'phase_sample', 'fitted_curve','z_grid','g1', 'fit_idx_1', 'fit_idx_2', 'corr_theory')

if 0
f = tight_subplot(1,3,[.06 .06],[.2 .1],[.1 .1]);

axes(f(1))
plot(z_grid, g1, 'o', 'MarkerSize', 4)
hold on
plot(z_grid, fitted_curve_1(z_grid), 'LineWidth',1.05);
plot(z_grid, corr_theory_1, 'LineWidth',1.05)
legend('Sampled (N = 10^4)', "Fitted (T = "+num2str(extracted_temperature_1*1e9, 3)+" nK)",...
    'Theory (Discrete)', 'FontSize', 10)
title('$T = 20 \; \rm nK$', 'Interpreter','latex')
legend box off
ylabel('$C(z)$', 'Interpreter', 'latex')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylim([0,1])

axes(f(2))
plot(z_grid, g2, 'o', 'MarkerSize',4)
hold on
plot(z_grid, fitted_curve_2(z_grid), 'LineWidth',1.05);
plot(z_grid, corr_theory_2, 'LineWidth',1.05)
legend('Sampled (N = 10^4)', "Fitted (T = "+num2str(extracted_temperature_2*1e9, 3)+ " nK)",...
    'Theory (Discrete)', 'FontSize', 10)
title('$T = 40 \; \rm nK$', 'Interpreter','latex')
legend box off
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylim([0,1])
yticks([])

axes(f(3))
plot(z_grid, g3, 'o',  'MarkerSize', 4)
hold on
plot(z_grid, fitted_curve_3(z_grid), 'LineWidth',1.05);
plot(z_grid, corr_theory_3, 'LineWidth',1.05)
legend('Sampled (N = 10^4)', "Fitted (T = "+num2str(extracted_temperature_3*1e9, 3)+ " nK)",...
    'Theory (Discrete)', 'FontSize', 10)
title('$T = 60\; \rm nK$', 'Interpreter','latex')
legend box off
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
ylim([0,1])
yticks([])

set(f, 'FontName', 'Times', 'FontSize', 14)
save('sample_thermometry_300microns', 'T1', 'T2', 'T3', 'phase_sample_1', 'phase_sample_2', 'phase_sample_3', 'g1', 'g2', 'g3', 'fitted_curve_1', 'fitted_curve_2', 'fitted_curve_3', 'z_grid',...)
'lambdaT1', 'lambdaT2', 'lambdaT3', 'hbar', 'm', 'kb','mean_density')
end
