clear all
close all
addpath('../../classes/')
addpath('../../plotting_func/')
addpath('results_shaking/scan7436/')
addpath('results_shaking/scan7440/')
addpath('results_shaking/scan7749/')
load('scan7436_exp_data_for_phase_extraction_all.mat')

%load density ripple data
density_ripple_all = data_structure_basic.density_profiles_full_stat;
mean_density = data_structure_basic.density_profiles;
evol_time = data_structure_basic.time;
z = data_structure_basic.z_axis;

n_max_fourier = 20;
deviation_mean_density = zeros(length(evol_time), length(z(20:65)));
fourier_mean_density_imag = zeros(length(evol_time), n_max_fourier);
fourier_mean_density_real = zeros(length(evol_time), n_max_fourier);
for j = 1:length(evol_time)
    dev = mean_density(20:65,j)- mean_density(20:65,1);
    fdev = fft(dev)';
    deviation_mean_density(j,:) = dev;
    fourier_mean_density_real(j,:) = real(fdev(2:n_max_fourier+1));
    fourier_mean_density_imag(j,:) = imag(fdev(2:n_max_fourier+1));

end

plot(evol_time*1e3, fourier_mean_density_imag(:,1).*1e-8, 'o-')
hold on
plot(evol_time*1e3, fourier_mean_density_imag(:,2).*1e-8, 'o-')
plot(evol_time*1e3, fourier_mean_density_imag(:,3).*1e-8, 'o-')
ylim([-1,1])
yline(0)
xline(0)

if 0
%load parameter

t_tof = 11.2e-3;


%set extraction parameter
n_max_fourier = 20;
idx_begin = 20; %The start of untruncated index
idx_end = 63; %The end of untruncated index
z_cut = z(idx_begin:idx_end);
shift = (z_cut(end) - z_cut(1))/2-z_cut(end);
z_cut = z_cut+shift;
condensate_length = z_cut(end)-z_cut(1);
dz = abs(z_cut(2)-z_cut(1));
mean_density_cut = mean_density(idx_begin:idx_end);
end