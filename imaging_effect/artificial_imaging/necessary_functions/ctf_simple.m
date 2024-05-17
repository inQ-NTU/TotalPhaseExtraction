function [ctf2D,defocus_arr_si,total_push_dist] = ctf_simple(extend,no_gridpoints,imaging_time,no_push_subdivision,number_photons,defocus_si,cloud_sigma_for_psf,numap)

%% constants

vrec = 5.8845e-3;

%% make a grid in nu-space

nu_spacing = 1/extend;
nu_grid_1D = get_fft_grid(nu_spacing,no_gridpoints);
[nuxg,nuyg] = meshgrid(nu_grid_1D,nu_grid_1D);
nutg = sqrt(nuxg.^2 + nuyg.^2);

%% calculate imaging push

timestep = imaging_time/no_push_subdivision; 
time_arr = timestep/2:timestep:imaging_time;

acc = number_photons*vrec/imaging_time; %acceleration due to absorbed photons
defocus_arr_si = acc*time_arr.^2/2;
total_push_dist = defocus_arr_si(end) - defocus_arr_si(1);

defocus_arr_si = defocus_arr_si - total_push_dist/2 - defocus_arr_si(1);
defocus_arr_si = defocus_arr_si + defocus_si;
        
%%

ctf2D = zeros(no_gridpoints,no_gridpoints,no_push_subdivision);
for push_ind=1:no_push_subdivision
    ctf2D(:,:,push_ind) = ifftshift(ctf_gauss_cloud(nutg,cloud_sigma_for_psf,defocus_arr_si(push_ind),numap));
end

end