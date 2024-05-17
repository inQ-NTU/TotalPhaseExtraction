function [atompics, backpics] = simulate_imaging( dens_2d, ...
                                                  grid_si, ...
                                                  binningsize, ...
                                                  avg_photon_cnts_per_pix, ...
                                                  ctf2D, ...
                                                  recoil_tf2D, ... 
                                                  recoil_flag, ...
                                                  shotnoise_flag, ...
                                                  camera_settings)
%%
% 
% input:
%   imaging_shift.....shift in gridpoints of the atoms, not in pixels
%   shift_dir.....can be 1 or 2, for Tandor 2 is the proper value


%% constants

h_si            = 6.6261e-34; %Placks constant
c               = 2.99e8; %speed of light
abs_cross       = 2.906e-13; % [m^2]
I_sat           = 16.6933; % Rb saturation intensity [W/m2]

E_ph            = h_si*c/780e-9;
alpha           = 1/camera_settings.alpha_fact;

%%  create artificial image

dx_grid             = ( grid_si(2) - grid_si(1) );
avgphotons          = avg_photon_cnts_per_pix * camera_settings.epc /(camera_settings.qe * binningsize^2); % avg photons per gridpoint
I0                  = avgphotons*E_ph/camera_settings.imaging_time/dx_grid^2; 

assert( abs(dx_grid*binningsize - camera_settings.pixelsize) < 1e-12)


no_push_subdivision = size(ctf2D, 1);
no_gridpoints       = length(grid_si);
I_atoms_PSF_arr     = zeros(no_push_subdivision, no_gridpoints, no_gridpoints);

if camera_settings.imaging_shift > 0
    shift_arr = linspace(0, camera_settings.imaging_shift/binningsize, no_push_subdivision);
    shift_arr = shift_arr -  mean(shift_arr);
    shift_arr = round(shift_arr);
end

% Exposure time is discretised in no_push_subdivision steps. For each step
% account for push of cloud and diffusion from photon re-emission  
for i = 1:no_push_subdivision
    
    if recoil_flag
        ctf_recoil_to_use = squeeze(recoil_tf2D(i,:,:));
        dens_2d_temp = ifft2(ctf_recoil_to_use.*fft2(dens_2d));
    else
        dens_2d_temp = dens_2d;
    end
    
    % calculate transmitted intensity for imaging pulse of plane wave
    I_atoms = I_sat/alpha .* lambertw(alpha*I0/I_sat .* exp(alpha*I0/I_sat - alpha*dens_2d_temp*abs_cross));
    
    % CTF must be applied to field (not intensity). As we are only
    % interested in the intensity, the complex part of the field does
    % not matter --> transform intensity to field via sqrt
    E_atoms         = sqrt(I_atoms);
    
    % apply PSF/CTF
    ctf2D_to_use    = squeeze(ctf2D(i,:,1));
    E_atoms_PSF     = ifft2(ctf2D_to_use.*fft2(E_atoms));
    
    % calculate intensity
    I_atoms_PSF_arr(i,:,:) = abs(E_atoms_PSF).^2;
    
    if camera_settings.imaging_shift > 0
        I_atoms_PSF_arr(i,:,:) = circshift(I_atoms_PSF_arr(i,:,:),shift_arr(i),camera_settings.shift_dir+1);
    end
end
 
I_atoms_avg     = squeeze(mean(I_atoms_PSF_arr,1));

% convert back to #photons
P_atoms         = I_atoms_avg*camera_settings.imaging_time/E_ph*dx_grid^2;
P_back          = ones(size(dens_2d)) * avgphotons; %background image = beam

% bin the image and account for quantum efficiency
P_atoms_binned  = binMatrix(P_atoms, binningsize) * camera_settings.qe;
P_back_binned   = binMatrix(P_back, binningsize) * camera_settings.qe;


if shotnoise_flag
    atompics   = double(poissrnd(P_atoms_binned))/camera_settings.epc;
    backpics   = double(poissrnd(P_back_binned))/camera_settings.epc;
else
    atompics   = P_atoms_binned/camera_settings.epc;
    backpics   = P_back_binned/camera_settings.epc;
end

% add background counts of camera
atompics = atompics + camera_settings.background_cnts;
backpics = backpics + camera_settings.background_cnts;

% rotate images to match typical output of experiment
atompics = permute(atompics, [2 1 3]);
backpics = permute(backpics, [2 1 3]);

end