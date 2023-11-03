clear all
close all
addpath('../classes/')
addpath('../plotting_func')

%This script is to give an example for single shot common phase and
%relative phase extraction

%Gas parameter
temp = 20e-9; %50 nK
n_fourier_cutoff = 15 ; %the maximum n 
pixnumz = 100;
t_tof = 15e-3;

flag_error_mitigation = 1;
flag_interaction_broadening = 0;
mean_density_profile = 'BoxPotential'; %Thomas-Fermi mean-density profile

if 0
%1. Initiate sampling with Bogoliubov sampling
sampling_suite = class_bogoliubov_phase_sampling(temp);
rel_phase = sampling_suite.generate_samples(n_fourier_cutoff, pixnumz);
com_phase = sampling_suite.generate_samples(n_fourier_cutoff, pixnumz);

%2. Fourier analysis of the input common phase
input_cosineCoeffs = sampling_suite.fourier_cosine_coeffs;
input_sineCoeffs = sampling_suite.fourier_sine_coeffs;
end

load('input_single_shot_phase.mat')
%%%%%3. Simulating TOF%%%%%%%
interference_suite = class_interference_pattern([rel_phase;com_phase], t_tof, mean_density_profile, flag_interaction_broadening);

x_grid = interference_suite.output_grid_x;
z_grid = interference_suite.output_grid_z;
insitu_density = interference_suite.insitu_density';
condensate_length = interference_suite.condensate_length_Lz;

rho_tof_full = interference_suite.tof_full_expansion();

%4. Compute density ripple
amp_full = trapz(x_grid, rho_tof_full, 2);

%calculating ripple function
%the ripple function is not defined in position where mean density is zero
%(near boundaries), so we introduce a cutoff to the boundary
boundary_factor = 0.05; %cut-down 5% of pixels from each side of the boundary
cut1 = floor(boundary_factor*pixnumz);
cut2 = floor((1-boundary_factor)*pixnumz);

ripple_func = 1 - amp_full(cut1:cut2)./(2*insitu_density(cut1:cut2));
%5. Use ripple information to deduce the spectrum of the common phase


%6. Error mitigation step
if flag_error_mitigation == 1
    %extract relative phase
    rel_suite = class_relative_phase_extraction(rho_tof_full, t_tof,flag_interaction_broadening);
    ext_rel_phase = rel_suite.fitting(rel_suite.init_phase_guess());

    %Compute density ripple without common phase
    interference_suite_woc = class_interference_pattern(ext_rel_phase, t_tof, mean_density_profile, flag_interaction_broadening); %woc -> without common phase
    rho_tof_full_woc = interference_suite_woc.tof_full_expansion();

    amp_woc = trapz(x_grid, rho_tof_full_woc, 2);
    ripple_error = 1 - amp_woc(cut1:cut2)./(2*insitu_density(cut1:cut2));
    ripple_func_new = ripple_func - ripple_error;
else
    ripple_func_new = ripple_func;
end

common_suite = class_common_phase_spectrum(ripple_func_new, z_grid(cut1:cut2), t_tof);
[output_cosineCoeffs, output_sineCoeffs] = common_suite.extract_com_spectrum(n_fourier_cutoff);
ext_com_phase = common_suite.extract_com_profile(z_grid);

%Plotting functions
z_grid = z_grid*1e6;
figure
g = tight_subplot(2,2,[.12 .12],[.1 .05],[.1 .05]);

axes(g(1))
plot(z_grid, amp_woc*1e-6, 'LineStyle','--', 'Color','black')
hold on
plot(z_grid, amp_full*1e-6,'Color','red')
plot(z_grid, 2*insitu_density*1e-6, 'Color','black')
xlabel('$z\; (\rm \mu m)$','Interpreter', 'latex')
ylabel('$n_{\rm TOF} \; (\rm \mu m^{-1})$', 'Interpreter', 'latex')
ax = gca;
ax.YAxis.Exponent = 2;
title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.95,0.8]);


axes(g(2))
plot(z_grid(cut1:cut2), ripple_func_new,'Color','red')
hold on
plot(z_grid(cut1:cut2), ripple_func,'black','LineStyle', '--')
ylabel('$\mathcal{R}(z,t)$', 'Interpreter', 'latex')
xlabel('$z\; \rm (\mu m)$', 'Interpreter', 'latex')
title('$\mathbf{b}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.95,0.8]);

axes(g(3))
k = 2*pi*(1:n_fourier_cutoff)/(condensate_length*1e6);
axes(g(3))
plot(k, input_cosineCoeffs,'Color','Black','LineStyle','-.')
hold on
plot(k, output_cosineCoeffs, 'o','Color','red')
ylb = ylabel('$a_k$','Interpreter','latex');
ylb.Position(1) = ylb.Position(1) + 0.02;
xlabel('$k\; (\rm \mu m^{-1})$', 'Interpreter', 'latex')
title('$\mathbf{c}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.95,0.8]);

f(1) = axes('Position',[.32 .15 .12 .12]);
box on
plot(k, input_sineCoeffs,'Color','Black','LineStyle','-.')
hold on
plot(k, output_sineCoeffs, 'o','MarkerSize',4,'Color','red')
ylb = ylabel('$b_k$','Interpreter', 'latex','FontSize', 16);
ylb.Position(1) = ylb.Position(1) + 0.1;
yticks([-0.2,0.2,0.6])

axes(g(4))
plot(z_grid, com_phase,'Color','Black','LineStyle', '-.')
hold on
plot(z_grid, ext_com_phase,'r')
ylabel('$\phi_+(z)$', 'Interpreter', 'latex')
xlabel('$z\; \rm (\mu m)$', 'Interpreter', 'latex')
set(g,'FontName','Times','FontSize',14)
set(f,'FontName','Times','FontSize',10)
title('$\mathbf{d}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.95,0.8]);
