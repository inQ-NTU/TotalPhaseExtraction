% Description:
% Example code for how to generate artificial images using the function
% create_artificial_images().
% The example will consist of 4 parts:
%   1. Generate some wavefunctions (here thermal equilibrium).
%   2. Simulate the corresponding 2D densities after TOF.
%   3. Simulate the imaging process.
%   4. Analyse the output images.

clear all; close all; clc;

% set path to all library functions
library_root = '../';
addpath(genpath(library_root))


%% Step 1: Generate/simulate wavefunctions to be imaged

% setup paramters
N_samples   = 5;            % number of shots/realizations 
N_atoms     = 8000;         % number of atoms
boxlength   = 100e-6;       % longitudinal boxtrap length
psf_DMD     = 1.0e-6;       % PSF of DMD responsible for boxtrap
T_si        = 50e-9;        % temperature
J_si        = 0;            % transversetunnel-coupling 
periodic_BC = false;        % flag for periodic boundary conditions
omegaT_si   = 2*pi*2e3;   % transver trapping frequency
N_wells     = 2;            % double well potential


% setup grid and calculate background density profile 
grid_dens     = linspace(-1.5*boxlength/2, 1.5*boxlength/2, 401);
dens_profile= density_TF_boxtrap(N_atoms/N_wells, boxlength, psf_DMD, grid_dens);

% sample thermal fluctuations and compute wavefunctions
psi         = generate_equilibrium_wavefunctions( N_samples, ...     
                                                  dens_profile, ...   
                                                  grid_dens, ...        
                                                  T_si, ...           
                                                  J_si, ...           
                                                  periodic_BC, ...    
                                                  omegaT_si, ...      
                                                  N_wells );  

%% Step 2: Calculate 2D density after TOF 

TOF_si          = 7e-3;
transverse_flag = false; %true; % if true, image via TAndor. Otherwise VAndor
fringe_spacing  = 15*1e-6;

[dens_2d, widths] = simulate_TOF( psi, ...
                                  TOF_si, ...
                                  grid_dens, ...
                                  fringe_spacing, ...
                                  omegaT_si, ...
                                  transverse_flag);


%% Step 3: Simulate artificial images

cloud_widths        = max(widths,[],1);     % width of each expanded cloud (needed for CTF)
imaging_intensity   = 0.25*16.6933;         % use 25% of saturation intensity
no_push_subdivisions= 20;                   % number of discretization steps in push simulation
imaging_system      = 'VAndor';             % string specifying the imaging system being simulated
shotnoise_flag      = true;                 % if true, account for photonic shotnoise
recoil_flag         = true;                 % if true, account for photon emmision recoil

[atompics, backpics, grid_pics] = create_artificial_images(  dens_2d, ...
                                                             grid_dens, ...
                                                             cloud_widths, ... 
                                                             imaging_intensity, ... 
                                                             no_push_subdivisions, ...
                                                             imaging_system, ... 
                                                             shotnoise_flag, ... 
                                                             recoil_flag );

%% Step 4: Analyse artificial images

% Create struct with abs_pic settings
abspic_settings.imaging_system  = imaging_system;
abspic_settings.background_correction = 1;  % 1 = back. correction on, 0 = off.
abspic_settings.useBackPics     = 2;
abspic_settings.useLoadedPics   = 0;
abspic_settings.gain            = 2;

abs_pic         = abs_pic_v2(abspic_settings);

density_imaged  = zeros(N_samples, length(grid_pics));

for i = 1:N_samples

    abs_pic.set_pics(atompics(:,:,i), backpics(:,:,i))
    abs_pic.load_data(); % creates absorption picture
    density_imaged(i,:) = abs_pic.get_density_y(); % get linescan
    
    figure
    abs_pic.plot_image()
end


% Compare density before and after imaging
density_raw = squeeze(trapz(grid_dens, dens_2d, 1))';

for i = 1:N_samples
    
    figure
    hold on
    box on
    
    plot(grid_dens*1e6, density_raw(i,:)*1e-6)
    plot(grid_pics*1e6, density_imaged(i,:)*1e-6)
    xlabel('z (um)')
    ylabel('density (atoms/um)')
    yline(0, '--')
    axis tight
    
    legend('Before imaging', 'After imaging')
    
    
end
