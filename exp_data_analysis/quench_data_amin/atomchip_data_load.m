%% Initialization
clear all; close all;


%% User input

scan_number = 8924; % Chose a scan number from the table in the slides.
t_ind = 1; % enter one of the time steps (from 1 to num_time_steps)
phi_lim_index = 1; % enter the index of the phi_lims = [pi, pi/2, pi/3, pi/4, pi/6]

%% Load data file

load(['data_filtered/scan_', int2str(scan_number), '_phase_data.mat'])

%% Description (please also see the slides or ask Amin) 

mean_density_per_um = saved.phase_mean_density; % MATLAB cell, averaged density for each time step in um^-1
phase_arr = saved.phase_arr; % MATLAB cell, all phases for each time step
density_arr = saved.density; % MATLAB cell, all densities for each time step
atom_number = saved.atom_number;
z_um = saved.z_um; % vector, grid points in um
z_l = saved.z_l; % int, index of the grid point for the left edge of the box
z_c = saved.z_c; % int, index of the grid point for the center of the box
z_r = saved.z_r; % int, index of the grid point for the right edge of the box
lT_um = saved.lambda_T_um; % Thermal coherence length in um
omega_perp_si = saved.omega_perp; % Transverse trapping frequency in SI units (s^-1, angular frequency)
times_ms = saved.times; % time vector in ms:
R2 = saved.R2;

% filtering 
phi_lims = saved.phi_lims; % phi_lims = [pi, pi/2, pi/3, pi/4, pi/6]
filtered_index = saved.phase_index_filtered;
filtered_out_index = saved.phase_index_filtered_out;
% [before_quench, right_after_the_quench, (time_evolution)]

num_time_steps = length(phase_arr); % number of time steps

num_shots_t_ind = size(phase_arr{t_ind}, 1); % number of shots for time step t_ind
% IMPORTANT NOTE: Don't be surprise if this number is lower than the count
% number in the scan overviec table (PowerPoint slides).The shots with atom
%numbers in the higher 15% or lower 15% of the distribution are discarded.

%% Shift the phase so that \phi{t}(:,z_c) is always between -pi to pi

phase_offset = repmat(angle( exp(1i*phase_arr{t_ind}(:,z_c)) )...
        - phase_arr{t_ind}(:,z_c)...
        ,1, size(phase_arr{t_ind},2));
    
shifted_phase_arr{t_ind} = phase_arr{t_ind} + phase_offset;

%%
clear saved

%% Print

fprintf('Scan %d has %d time steps.\n', scan_number, num_time_steps)
fprintf('Scan %d - time step %d has %d shots.\n', scan_number, ...
                                                    t_ind, num_shots_t_ind)
%% Example plots

% plot atom number for all shots time index t_ind
figure('Name','Atom number')
clf

plot(atom_number{t_ind}, '-o')
axis([-inf inf 0 inf])
title({['Scan', int2str(scan_number)], ...
    ['Atom number, t = ' num2str(times_ms(t_ind), 3), ' ms']})
ylabel('N');
xlabel('shot number');

% plot mean density for time index t_ind
figure('Name','Averaged density')
clf

plot(z_um, mean_density_per_um{t_ind})
axis([-inf inf -inf inf])
xline(z_um(z_l)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_c)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_r)); % Comment out if using MATLAB older than R2018b!
title({['Scan', int2str(scan_number)], ...
    ['Averaged density, t = ' num2str(times_ms(t_ind), 3), ' ms']})
ylabel('\rho(z) (\mu m^{-1})');
xlabel('z (\mu m)');

% plot all density profiles for the time index t_ind
figure('Name','All density profiles')
clf
hold on

for ii = 1:size(phase_arr{t_ind},1)
    plot(z_um, squeeze(density_arr{t_ind}(ii,:)))
end

axis([-inf inf 0 inf])
xline(z_um(z_l)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_c)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_r)); % Comment out if using MATLAB older than R2018b!
title({['Scan', int2str(scan_number)], ...
    ['All density profiles, t = ' num2str(times_ms(t_ind), 3), ' ms']})
ylabel('\rho(z) (1/ \mu m)');
xlabel('z (\mu m)');


% plot all phase profiles for the time index t_ind
figure('Name','All phase profiles')
clf
hold on

for ii = 1:size(phase_arr{t_ind},1)
    plot(z_um, squeeze(phase_arr{t_ind}(ii,:)))
end

axis([-inf inf -inf inf])
xline(z_um(z_l)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_c)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_r)); % Comment out if using MATLAB older than R2018b!
title({['Scan', int2str(scan_number)], ...
    ['All phase profiles, t = ' num2str(times_ms(t_ind), 3), ' ms']})
ylabel('\phi(z)');
xlabel('z (\mu m)');

%%
figure('Name','All shifted phase profiles')
clf
hold on

for ii = 1:size(phase_arr{t_ind},1)
    plot(z_um, squeeze(shifted_phase_arr{t_ind}(ii,:)))
end

axis([-inf inf -inf inf])
xline(z_um(z_l)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_c)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_r)); % Comment out if using MATLAB older than R2018b!
title({['Scan', int2str(scan_number)], ...
    ['All shifted phase profiles, t = ' num2str(times_ms(t_ind), 3), ' ms']})
ylabel('\phi(z)');
xlabel('z (\mu m)');
%%

% plot filtered OUT phase profiles
if size(filtered_out_index{t_ind, phi_lim_index},1) > 1
    figure('Name','Filtered OUT phase profiles')
    clf
    hold on

    for ii = 1:size(filtered_out_index{t_ind, phi_lim_index},1)
        plot(z_um, squeeze(phase_arr{t_ind}(filtered_out_index{t_ind, phi_lim_index},:)))
    end

    axis([-inf inf -inf inf])
    xline(z_um(z_l)); % Comment out if using MATLAB older than R2018b!
    xline(z_um(z_c)); % Comment out if using MATLAB older than R2018b!
    xline(z_um(z_r)); % Comment out if using MATLAB older than R2018b!
    title({['Scan', int2str(scan_number)], ...
        ['Filtered out phase profiles, \phi_{lim} = ' num2str(phi_lims(phi_lim_index), 3) ', t = ' num2str(times_ms(t_ind), 3), ' ms']})
    ylabel('\phi(z)');
    xlabel('z (\mu m)');

end
% plot all R^2 for the time index t_ind
figure('Name','All R2')
clf
hold on

for ii = 1:size(R2{t_ind},1)
    plot(z_um, squeeze(R2{t_ind}(ii,:)))
end

axis([-inf inf 0 1])
xline(z_um(z_l)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_c)); % Comment out if using MATLAB older than R2018b!
xline(z_um(z_r)); % Comment out if using MATLAB older than R2018b!
title({['Scan', int2str(scan_number)], ...
    ['All R^2, t = ' num2str(times_ms(t_ind), 3), ' ms']})
ylabel('R^2');
xlabel('z (\mu m)');
% end