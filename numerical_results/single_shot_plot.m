close all;
clear all;

addpath('../classes/')
addpath('../plotting_func/')

mean_density = 75e6;
temperature_com = 20e-9;
temperature_rel = 30e-9;
scaling_factor = sqrt(2);
max_n_fourier = 20;
pixnum = 201;
t_tof = 11e-3;

load('single_shot_input.mat')

if 0
com_suite = class_bogoliubov_sampling(temperature_com, mean_density, scaling_factor);
rel_suite = class_bogoliubov_sampling(temperature_rel, mean_density, scaling_factor);

[com_phase, com_density_fluct] = com_suite.generate_fluct_samples(max_n_fourier, pixnum);
[rel_phase, rel_density_fluct] = rel_suite.generate_fluct_samples(max_n_fourier, pixnum);
end

phase_1 = (com_phase + rel_phase)/2;
phase_2 = (com_phase - rel_phase)/2;

delta_rho_1 = (com_density_fluct + rel_density_fluct)/2;
delta_rho_2 = (com_density_fluct - rel_density_fluct)/2;

insitu_density_1 = mean_density*ones(1,pixnum)+delta_rho_1;
insitu_density_2 = mean_density*ones(1,pixnum)+delta_rho_2;

interference_suite = class_interference_pattern([phase_1;phase_2], [insitu_density_1; insitu_density_2], t_tof);

rho_tof = interference_suite.tof_full_expansion();
x_axis  = interference_suite.output_grid_x;
z_axis = interference_suite.output_grid_z;
y_axis = x_axis;
sigma_yx = 30e-6;

density_ripple = trapz(x_axis, rho_tof, 2);
density_ripple_2d = density_ripple*exp(-y_axis.^2/(sigma_yx)^2);

bulk_start = -40e-6;
bulk_end = 40e-6;

[~, idx_start] = min(abs(z_axis - bulk_start));
[~, idx_end] = min(abs(z_axis - bulk_end));

cut_z_axis = z_axis(idx_start:idx_end);
cut_density_ripple = density_ripple(idx_start:idx_end);

ext_com_suite = class_common_phase_extraction(cut_density_ripple, 2*mean_density*ones(1, idx_end - idx_start+1), ...
    cut_z_axis, t_tof);

[cosineCoeffs, sineCoeffs] = ext_com_suite.extract_com_spectrum(max_n_fourier);
ext_com_phase = ext_com_suite.reconstruct_com_phase(cosineCoeffs, sineCoeffs);


figure
imagesc(rho_tof);
colormap(gge_colormap)
view([180,90,90])
xticks([])
yticks([])
axis off
grid off

figure
imagesc(density_ripple_2d')
colormap(gge_colormap)
xticks([])
yticks([])
axis off
grid off

figure
f = tight_subplot(2,1, [0.1, 0.1], [0.15, 0.1], [0.15,0.15]);
axes(f(1))
plot(z_axis*1e6, density_ripple*1e-6, 'Color', 'red', 'LineWidth',1.05)
xticks([])
yline(2*mean_density*1e-6, '-.', 'LineWidth',1.05)
ylabel('$n_{\rm tof}\; \rm (\mu m^{-1})$', 'Interpreter','latex')
xline(bulk_start*1e6, '--', 'LineWidth',1.05)
xline(bulk_end*1e6, '--', 'LineWidth',1.05)
ax = gca;
ax.YAxis.Exponent = 2;
ax.LineWidth = 1.1;

axes(f(2))
plot(z_axis*1e6, com_phase, '-.', 'Color','black', 'LineWidth',1.05)
hold on
plot(cut_z_axis*1e6, ext_com_phase, 'Color','red', 'LineWidth',1.05)
ylabel('$\phi_+(z)$', 'Interpreter','latex')
xlabel('$z\; \rm (\mu m)$', 'Interpreter','latex')
xline(bulk_start*1e6, '--', 'LineWidth',1.05)
xline(bulk_end*1e6, '--', 'LineWidth',1.05)
ax = gca; 
ax.LineWidth = 1.1;

set(f, 'FontName', 'Times', 'FontSize', 20)

%save('single_shot_input.mat', 'com_phase', 'rel_phase', 'com_density_fluct', 'rel_density_fluct', 'mean_density', 't_tof', ...
%    'rho_tof', 'density_ripple', 'density_ripple_2d', 'idx_end','idx_start', 'bulk_end', 'bulk_start', 'ext_com_phase')