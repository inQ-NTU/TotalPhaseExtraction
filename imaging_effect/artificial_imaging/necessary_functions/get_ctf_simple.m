function [ctf2D,defocus_arr_si] = get_ctf_simple(input_data)

% user input
grid_si = input_data.grid_si;
binningsize = input_data.binningsize;
avg_photon_cnts_per_pix = input_data.avg_photon_cnts_per_pix;
epc = input_data.epc;
qe = input_data.qe;
alpha_fact = input_data.alpha_fact;
no_push_subdivision = input_data.no_push_subdivision;
push_flag = input_data.push_flag;
imaging_time = input_data.imaging_time;
cloud_sigma_for_psf = input_data.cloud_sigma_for_psf;
defocus_um = input_data.defocus_um;
numap = input_data.numap;

%% calculate some helpful quantities

no_gridpoints = length(grid_si);
grid_spacing_si = grid_si(2) - grid_si(1);
extend = no_gridpoints * grid_spacing_si;

%% calculate number of scattered photons per atom

number_photons = get_no_ph(grid_spacing_si,binningsize,avg_photon_cnts_per_pix,epc,qe,alpha_fact);

%% calculate

if push_flag
    [ctf2D,defocus_arr_si] = ctf_simple(extend,no_gridpoints,imaging_time,no_push_subdivision,number_photons,defocus_um*1e-6,cloud_sigma_for_psf,numap);
else
    no_push_subdivision_temp = 1;
    [ctf2D,defocus_arr_si] = ctf_simple(extend,no_gridpoints,imaging_time,no_push_subdivision_temp,number_photons,defocus_um*1e-6,cloud_sigma_for_psf,numap);
    ctf2D = repmat(ctf2D,1,1,no_push_subdivision);
end

end