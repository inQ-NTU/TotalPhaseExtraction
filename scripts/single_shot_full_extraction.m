clear all
close all
addpath('../classes/')
addpath('../plotting_func')

%This script is to give an example for single shot common phase and
%relative phase extraction

%gas parameter
temp = 50e-9; %50 nK
n_fourier_cutoff = 15; %the maximum n 
pixel_size = 1e-6; %1 microns 
condensate_length = 100e-6; 
pixnumz = floor(condensate_length/pixel_size); %number of longitudinal pixels
z_grid = linspace(-condensate_length/2, condensate_length/2-pixel_size, pixnumz);
x_grid = z_grid;

flag_interaction_broadening = 1; %turn-on/off interaction broadening
mean_density_profile = 'InverseParabola'; %Thomas-Fermi mean-density profile
t_tof = 15e-3;

%1. Initiate sampling with Ornstein-Uhlenbeck (OU) process
%Cite references
%if 0
OU_suite = class_OU_phase_sampling(temp);
rel_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz);
%rel_phase = zeros(1,pixnumz);
%com_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz);
com_phase = zeros(1,pixnumz); 
%com_phase = pi*cos(6*pi*z_grid/condensate_length);
input_cosineCoeffs = OU_suite.fourier_cosine_coeffs;
input_sineCoeffs = OU_suite.fourier_sine_coeffs;
%end
%load('single_shot_phase_profile.mat')
%%%%%Simulating TOF%%%%%%%
interference_suite = class_interference_pattern([rel_phase;com_phase],t_tof, mean_density_profile, flag_interaction_broadening);

rho_tof_trans = interference_suite.tof_transversal_expansion();
rho_tof_full = interference_suite.tof_full_expansion();

%Normalize to 10,000 atoms
rho_tof_trans = interference_suite.normalize(rho_tof_trans, 10^4);
rho_tof_full = interference_suite.normalize(rho_tof_full, 10^4);

%Compute density ripple
amp_trans = trapz(x_grid, rho_tof_trans, 2);
amp_full = trapz(x_grid, rho_tof_full, 2);

%calculating ripple function
%the ripple function is not defined in position where mean density is zero
%(near boundaries), so we introduce a cutoff to the boundary
boundary_factor = 0.1; %cut-down 10% of pixels from each side of the boundary
cut1 = boundary_factor*pixnumz;
cut2 = (1-boundary_factor)*pixnumz;

ripple_func = 1 - amp_full(cut1:cut2)./amp_trans(cut1:cut2);

%Use ripple information to deduce the spectrum of the common phase
common_suite = class_common_phase_spectrum(ripple_func, z_grid(cut1:cut2), t_tof);
[output_cosineCoeffs, output_sineCoeffs] = common_suite.extract_com_spectrum(n_fourier_cutoff);
ext_com_phase = common_suite.extract_com_profile(z_grid);


%Plotting functions
figure
g = tight_subplot(2,2,[.15 .12],[.15 .05],[.1 .05]);
axes(g(1))
imagesc(z_grid*1e6, x_grid*1e6, rho_tof_full')
xlabel('$z\; (\rm \mu m)$','Interpreter', 'latex')
ylabel('$x\; (\rm \mu m)$','Interpreter', 'latex')
%ylb = ylabel('$x \; (\rm \mu m)$', 'Interpreter','latex');
%ylb.Position(1) = ylb.Position(1) - 0.05;
%title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.95,0.8]);
colormap(gge_colormap)
yticks([-50,0,50])
xticks([-50,0,40])
xticklabels([-50,0,50])

axes(g(2))
plot(z_grid*1e6, amp_full*1e-6,'Color','black')
hold on
plot(z_grid*1e6, amp_trans*1e-6, '--','Color','black')
xlabel('$z\; (\rm \mu m)$','Interpreter', 'latex')
ylabel('$n_{\rm TOF} \; (\rm \mu m^{-1})$', 'Interpreter', 'latex')
ax = gca;
ax.YAxis.Exponent = 2;

axes(g(3))
k = 2*pi*(1:n_fourier_cutoff)/(condensate_length*1e6);
axes(g(3))
%plot(k, input_cosineCoeffs,'Color','Black')
%hold on
plot(k, output_cosineCoeffs, 'o','Color','red')
ylb = ylabel('$a_k$','Interpreter','latex');
ylb.Position(1) = ylb.Position(1) + 0.02;
xlabel('$k\; (\rm \mu m^{-1})$', 'Interpreter', 'latex')
ax = gca;
ax.YAxis.Exponent = -1;

%f(1) = axes('Position',[.28 .34 .1 .1]);
%box on
%plot(k, input_sineCoeffs,'Color','Black')
%hold on
%plot(k, output_sineCoeffs, 'o','MarkerSize',4,'Color','red')
%ylb = ylabel('$b_k$','Interpreter', 'latex','FontSize', 16);
%ylb.Position(1) = ylb.Position(1) + 0.1;
%yticks([0,0.4])
%title('$\mathbf{a}$','FontName','Times','Color','black','Units', 'normalized','Interpreter','latex','Position',[0.95,0.8]);

axes(g(4))
plot(z_grid*1e6, com_phase, '--','Color','Black')
hold on
plot(z_grid*1e6, ext_com_phase,'Color','red')
ylabel('$\phi_+(z)$', 'Interpreter', 'latex')
xlabel('$z\; \rm (\mu m)$', 'Interpreter', 'latex')
%yticks([-1,0,1])
set(g,'FontName','Times','FontSize',16)