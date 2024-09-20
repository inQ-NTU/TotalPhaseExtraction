clear all; close all;
addpath('../../../classes/')
addpath('../../../plotting_func/')
addpath('../results_shaking/scan7436/')
load('scan7436_exp_data_for_phase_extraction_all.mat')

%load density ripple data
density_ripple_all = data_structure_basic.density_profiles_full_stat;
mean_density = data_structure_basic.density_profiles(:,1);

%load parameter
evol_time = data_structure_basic.time;
z = data_structure_basic.z_axis;

%set extraction parameter
%n_max_fourier = 20;
idx_begin = 20; %The start of untruncated index
idx_end = 64; %The end of untruncated index
z_cut = z(idx_begin:idx_end);
shift = (z_cut(end) - z_cut(1))/2-z_cut(end);
z = (z+shift)*1e6;
condensate_length = data_structure_basic.L_fit;
mean_density_cut = mean_density(idx_begin:idx_end);

thermal_idx = 1;
shaken_idx_1 = 7;
shaken_idx_2 = 9;
shaken_idx_3 = 11;
shaken_idx_4 = 13; 
shaken_idx_5 = 15;

for j = 1:length(evol_time)
    density_ripple_thermal = 2*density_ripple_all{thermal_idx}.*1e-6;
    density_ripple_shaken_1 = 2*density_ripple_all{shaken_idx_1}.*1e-6;
    density_ripple_shaken_2 = 2*density_ripple_all{shaken_idx_2}.*1e-6;
    density_ripple_shaken_3 = 2*density_ripple_all{shaken_idx_3}.*1e-6;
    density_ripple_shaken_4 = 2*density_ripple_all{shaken_idx_4}.*1e-6;
    density_ripple_shaken_5 = 2*density_ripple_all{shaken_idx_5}.*1e-6;
end

mean_density_ripple_thermal = mean(density_ripple_thermal,2);
mean_density_ripple_shaken_1 = mean(density_ripple_shaken_1,2);
mean_density_ripple_shaken_2 = mean(density_ripple_shaken_2,2);
mean_density_ripple_shaken_3 = mean(density_ripple_shaken_3,2);
mean_density_ripple_shaken_4 = mean(density_ripple_shaken_4,2);
mean_density_ripple_shaken_5 = mean(density_ripple_shaken_5,2);

mean_density = mean(mean_density_ripple_thermal(idx_begin:idx_end));
h = tight_subplot(3,2, [0.1, 0.08], [0.1, 0.05], [0.12, 0.05]);
axes(h(1))
plot(z,density_ripple_thermal, 'Color',[0.7,0.7,0.7])
hold on
plot(z,mean_density_ripple_thermal, 'Color','red','LineWidth',1.1)
xline(z(idx_begin), '--')
xline(z(idx_end), '--')
yline(mean_density)
xlim([-45,45])
xticks([])
ylim([0,250])
ylabel('$n_{\rm tof}(\rm \mu m^{-1})$', 'Interpreter','latex')
title('$t = -30\; \rm ms$', 'Interpreter','latex')

axes(h(2))
plot(z,density_ripple_shaken_1, 'Color',[0.7,0.7,0.7])
hold on
plot(z,mean_density_ripple_shaken_1, 'Color','red', 'LineWidth',1.1)
xline(z(idx_begin), '--')
xline(z(idx_end), '--')
xlim([-45,45])
xticks([])
ylim([0,250])
yticks([])
yline(mean_density)
title('$t = 0\; \rm ms$', 'Interpreter','latex')

axes(h(3))
plot(z,density_ripple_shaken_2, 'Color',[0.7,0.7,0.7])
hold on
plot(z,mean_density_ripple_shaken_2, 'Color','red', 'LineWidth',1.1)
xline(z(idx_begin), '--')
xline(z(idx_end), '--')
xlim([-45,45])
xticks([])
ylim([0,250])
ylabel('$n_{\rm tof}(\rm \mu m^{-1})$', 'Interpreter','latex')
yline(mean_density)
title('$t = 10\; \rm ms$', 'Interpreter','latex')

axes(h(4))
plot(z,density_ripple_shaken_3, 'Color',[0.7,0.7,0.7])
hold on
plot(z,mean_density_ripple_shaken_3, 'Color','red', 'LineWidth',1.1)
xline(z(idx_begin), '--')
xline(z(idx_end) ,'--')
xlim([-45,45])
xticks([])
ylim([0,250])
yticks([])
yline(mean_density)
title('$t = 20\; \rm ms$', 'Interpreter','latex')

axes(h(5))
plot(z,density_ripple_shaken_4, 'Color',[0.7,0.7,0.7])
hold on
plot(z,mean_density_ripple_shaken_4, 'Color','red', 'LineWidth',1.1)
xline(z(idx_begin), '--')
xline(z(idx_end), '--')
xlim([-45,45])
ylim([0,250])
ylabel('$n_{\rm tof}(\rm \mu m^{-1})$', 'Interpreter','latex')
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
yline(mean_density)
title('$t = 30\; \rm ms$', 'Interpreter','latex')

axes(h(6))
plot(z,density_ripple_shaken_5, 'Color',[0.7,0.7,0.7])
hold on
plot(z,mean_density_ripple_shaken_5, 'Color','red', 'LineWidth',1.1)
xline(z(idx_begin), '--')
xline(z(idx_end), '--')
xlim([-45,45])
ylim([0,250])
yticks([])
xlabel('$z\; (\rm \mu m)$', 'Interpreter','latex')
yline(mean_density)
title('$t = 40\; \rm ms$', 'Interpreter','latex')

set(h, 'FontName', 'Times', 'FontSize', 14)