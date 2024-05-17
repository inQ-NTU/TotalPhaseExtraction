function [atompics, backpics, grid_pics] = create_artificial_images(  dens_2d, ... % 2d density profiles on square grid
                                                                      grid_dens, ... % underlying grid (1d)
                                                                      cloud_widths, ... % widths of clouds (if all clouds have same width, pass double)
                                                                      imaging_intensity, ... % average intensity of imaging light in SI-units
                                                                      no_push_subdivisions, ... % number of discretization steps in push simulation
                                                                      imaging_system, ... % string specifying the imaging system being simulated
                                                                      shotnoise_flag, ... % if true, account for photonic shotnoise
                                                                      recoil_flag ) % if true, account for photon emmision recoil
% ---------------------------------------------------------------
%
% Creates artificial images from simulated profiles.
%
% Original code by Thomas Schweigler 2016-04
% Update by Frederik Moller 2019-03
%   - now takes saturation into account
%   - works for both VAndor and TAndor
%   - works for both single- and double-well
% Update by Frederik Moller 2022-02
%   - moved image plotting and saving to seperate functions
%   - made inputs more intuitive


%% Set camera properties

if strcmpi(imaging_system,'verticalAndor')|| strcmpi(imaging_system,'VAndor')
    camera_settings.pixelsize               = 1.9465e-06;   % pixel size in SI-units
    camera_settings.imaging_time            = 50e-6;        % exposure time in SI-units
    camera_settings.avg_cnts_per_pix_total  = 2300;         % (including background counts)
    camera_settings.epc                     = 1.84;         % electrons per count
    camera_settings.qe                      = 0.8;          % quantum efficiency
    camera_settings.background_cnts         = 1001;         % static background counts of camera
    camera_settings.alpha_fact              = 1;            % reduction of scattering cross-section from imaging light polarization
    camera_settings.numap                   = 0.079;        % numerical aperture
    camera_settings.push_flag               = false;        % do not account for push in VAndor images
    camera_settings.defocus_um              = 32.7;          % defocus in micrometers
    camera_settings.imaging_shift           = 0;            % number of pixels that shift on CCD during exposure
    camera_settings.shift_dir               = 1;            % direction that pixels shift in
    
elseif strcmpi(imaging_system,'transversalAndor')|| strcmpi(imaging_system,'TAndor')
    camera_settings.pixelsize               = 1.0492e-6;    % pixel size in SI-units
    camera_settings.imaging_time            = 75e-6;        % exposure time in SI-units
    camera_settings.avg_cnts_per_pix_total  = 765;          % (including background counts)
    camera_settings.epc                     = 2.87;         % electrons per count
    camera_settings.qe                      = 0.57;         % quantum efficiency
    camera_settings.background_cnts         = 458;          % static background counts of camera
    camera_settings.alpha_fact              = 1/0.54;       % reduction of scattering cross-section from imaging light polarization
    camera_settings.numap                   = 0.18;         % numerical aperture
    camera_settings.push_flag               = true;         % account for push in TAndor images
    camera_settings.defocus_um              = -20;          % defocus in micrometers
    camera_settings.imaging_shift           = 0;            % NEED TO FIX THIS!
    camera_settings.shift_dir               = 1;            % NEED TO FIX THIS!
    
else
    error('Unrecognised imaging system!')
end
%% Setup grids for interpolating 2d density to image grid

% Create simulation grid comensurate with number of pixels (for better results
% use a finer grid for imagining simulation and then bin to camera pixels)

L_img       = grid_dens(end)-grid_dens(1); % side length of image in meters
N_pix       = ceil(L_img/camera_settings.pixelsize); % number of pixels needed to cover grid
binningsize = ceil(length(grid_dens)/N_pix); % number of gridpoints per pixel

grid_simul  = (1:N_pix*binningsize)*camera_settings.pixelsize/binningsize;
grid_pics   = camera_settings.pixelsize*(1:N_pix);

% make sure all grids are centered on 0
grid_dens   = grid_dens - mean(grid_dens); 
grid_simul  = grid_simul - mean(grid_simul);
grid_pics   = grid_pics - mean(grid_pics);

%% Calculate ctf (coherent transfer function)

% Calculate average number of photon counts per camera pixel. If no
% imaging intensity is given (it is zero or less) use default value.
if imaging_intensity > 0
    h_si                    = 6.6261e-34;   % Placks constant
    c                       = 2.99e8;       % speed of light
    E_ph                    = h_si*c/780e-9;% photon energy

    px_area                 = camera_settings.pixelsize^2; % area of one pixel in SI-units (at the plane of the atoms)
    photon_density          = imaging_intensity*camera_settings.imaging_time/E_ph; % total photon density in imaging light
    avg_photon_cnts_per_pix = photon_density*px_area*(camera_settings.qe/camera_settings.epc);
else
    % default value based on measurement
    avg_photon_cnts_per_pix = camera_settings.avg_cnts_per_pix_total - camera_settings.background_cnts;
end

% Pack parameters in struct
ctf_params = camera_settings;
ctf_params.binningsize = binningsize;
ctf_params.no_push_subdivision = no_push_subdivisions;
ctf_params.grid_si = grid_simul;
ctf_params.avg_photon_cnts_per_pix = avg_photon_cnts_per_pix;

assert(isvector(cloud_widths), 'Only pass single cloud width for each realisation.')

ctf2D = zeros(  no_push_subdivisions, ... 
                length(grid_simul), ...
                length(grid_simul), ...
                length(cloud_widths));
recoil_tf2D = ctf2D;

% Calculate coherent transfer functions
for i = 1:length(cloud_widths)
    ctf_params.cloud_sigma_for_psf = cloud_widths(i);

    ctf2D(:,:,:,i) = permute(get_ctf_simple(ctf_params), [3 1 2]);
    
    if recoil_flag
        recoil_tf2D(:,:,:,i) = get_tf_recoil(ctf_params);
    end
end


%% create artificial images

N_pics      = size(dens_2d, 3);
atompics    = zeros(N_pix, N_pix, N_pics );
backpics    = zeros(N_pix, N_pix, N_pics );

parfor i = 1:N_pics       
    % interpolate to pixel-commensurate grid
    dens_atoms = interp2(grid_dens', grid_dens, dens_2d(:,:,i), grid_simul', grid_simul, 'cubic', 0 ); % all extrapolated values are 0
    
    % Calculate the number of counts distributed over the camera chip
    [P_atoms, P_noatoms] = simulate_imaging( dens_atoms, ...
                                             grid_simul, ...
                                             binningsize, ...
                                             avg_photon_cnts_per_pix, ...
                                             ctf2D(:,:,:,i), ...
                                             recoil_tf2D(:,:,:,i), ... 
                                             recoil_flag, ...
                                             shotnoise_flag, ...
                                             camera_settings);

                                                    
    atompics(:,:,i) = P_atoms;
    backpics(:,:,i) = P_noatoms;

end

end %ArtificialImaging