close all
clear all
addpath('../classes')
%Check sample average
%gas parameter
temp =50e-9; %50 nK
n_fourier_cutoff = 15; %the maximum n 
pixel_size = 1e-6; %1 microns 
condensate_length = 100e-6; 
pixnumz = floor(condensate_length/pixel_size); %number of longitudinal pixels
z_grid = linspace(-condensate_length/2, condensate_length/2-pixel_size, pixnumz);
x_grid = z_grid;
N_samples = 1000;
t_tof = 15e-3;
flag_interaction_broadening = 1;

%1. Initiate sampling with Ornstein-Uhlenbeck (OU) process
%Cite references
OU_suite = class_OU_phase_sampling(temp);
rel_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz, N_samples);
com_phase = OU_suite.generate_OU_samples(n_fourier_cutoff, pixnumz, N_samples);

com_phase_input_cosineCoeffs = OU_suite.fourier_cosine_coeffs;
com_phase_input_sineCoeffs = OU_suite.fourier_sine_coeffs;
both_input = com_phase_input_cosineCoeffs + 1j*com_phase_input_sineCoeffs;

mean_cosineCoeffs_input = zeros(1, n_fourier_cutoff);
mean_sineCoeffs_input = zeros(1,n_fourier_cutoff);
mean_both_input = zeros(1,n_fourier_cutoff);
for n = 1:n_fourier_cutoff
    mean_cosineCoeffs_input(n) = mean(com_phase_input_cosineCoeffs(:,n).^2);
    mean_sineCoeffs_input(n) = mean(com_phase_input_sineCoeffs(:,n).^2);
    mean_both_input(n) = mean(abs(both_input(:,n)).^2);
end

com_phase_output_cosineCoeffs = zeros(N_samples, n_fourier_cutoff); 
com_phase_output_sineCoeffs = zeros(N_samples, n_fourier_cutoff);
output_both = zeros(N_samples, n_fourier_cutoff);
count = 0;
for i = 1:N_samples
    %%%%%Simulating TOF%%%%%%%
    interference_suite = class_interference_pattern([rel_phase(i,:);com_phase(i,:)],t_tof, 'BoxPotential', flag_interaction_broadening);

    rho_tof_trans = interference_suite.tof_transversal_expansion();
    rho_tof_full = interference_suite.tof_full_expansion();

    %Normalize to 10,000 atoms
    %rho_tof_trans = interference_suite.normalize(rho_tof_trans, 10^4);
    %rho_tof_full = interference_suite.normalize(rho_tof_full, 10^4);

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
    com_phase_output_cosineCoeffs(i,:) = output_cosineCoeffs;
    com_phase_output_sineCoeffs(i,:) = output_sineCoeffs;
    output_both(i,:) = output_cosineCoeffs + 1j*output_sineCoeffs;
    count = count+1;
    disp(count)
end

mean_cosineCoeffs_output = zeros(1, n_fourier_cutoff);
mean_sineCoeffs_output = zeros(1,n_fourier_cutoff);
mean_magnitude = zeros(1,n_fourier_cutoff);


for n = 1:n_fourier_cutoff
    mean_cosineCoeffs_output(n) = mean(com_phase_output_cosineCoeffs(:,n).^2);
    mean_sineCoeffs_output(n) = mean(com_phase_output_sineCoeffs(:,n).^2);
    mean_magnitude(n) = mean(abs(output_both(:,n)).^2);
end

%k = (2*pi*(2:15)/(condensate_length*1e6));
n_list = (2:15);
% Define Start points, fit-function and fit curve
fitfun = fittype( @(a,b,x) a./(b+x.^2));
fitted_curve_out = fit(n_list',mean_magnitude(2:15)',fitfun);
fitted_curve_in = fit(n_list', mean_both_input(2:15)', fitfun);
% Save the coeffiecient values for a,b,c and d in a vector
coeffvals1 = coeffvalues(fitted_curve_out)
coeffvals2 = coeffvalues(fitted_curve_in)
plot(n_list, mean_both_input(2:end),'x')
hold on
plot(n_list, mean_magnitude(2:end),'o')
plot(n_list, fitted_curve_out(n_list))

m = OU_suite.m;
kb = OU_suite.kb;
n0 = OU_suite.max_longitudinal_density;
h = OU_suite.hbar*2*pi;

temp1 = n0*h^2*coeffvals1(1)/(2*m*kb*condensate_length)
temp2 = n0*h^2*coeffvals2(1)/(2*m*kb*condensate_length)

disp(abs(temp1 - temp2)/temp2)
