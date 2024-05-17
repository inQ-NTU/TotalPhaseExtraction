function absorption_image_dens = absorption_imaging(dens_2d, grid_dens, cloud_width)
    %Imaging setup
    imaging_intensity   = 0.25*16.6933;         % use 25% of saturation intensity
    no_push_subdivisions= 20;                   % number of discretization steps in push simulation
    imaging_system      = 'TAndor';             % string specifying the imaging system being simulated
    shotnoise_flag      = true;                 % if true, account for photonic shotnoise
    recoil_flag         = true;                 % if true, account for photon emmision recoil

    [img_dens, backpics] = create_artificial_images(dens_2d, ...
                                 grid_dens, ...
                                 cloud_width, ... 
                                 imaging_intensity, ... 
                                 no_push_subdivisions, ...
                                 imaging_system, ... 
                                 shotnoise_flag, ... 
                                 recoil_flag);


    %%Analyse artificial images
    % Create struct with abs_pic settings
    abspic_settings.imaging_system  = imaging_system;
    abspic_settings.background_correction = 1;  % 1 = back. correction on, 0 = off.
    abspic_settings.useBackPics     = 2;
    abspic_settings.useLoadedPics   = 0;
    abspic_settings.gain            = 2;

    abs_pic         = abs_pic_v2(abspic_settings);
    abs_pic.set_pics(img_dens(:,:), backpics(:,:))
    abs_pic.load_data(); % creates absorption picture
    absorption_image_dens = abs_pic.gain.*256.*abs_pic.ROI_data;
    %absorption_image_dens = img_dens;
end
