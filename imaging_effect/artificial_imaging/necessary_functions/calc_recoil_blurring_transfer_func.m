function [transfer_func,dist_2D_conv_norm] = calc_recoil_blurring_transfer_func(time,avg_photon_cnts_per_pix,epc,qe,rg,grid_spacing,binningsize,alpha_fact)

%% constants

absCross=2.906e-13; %m^2

%% calculate number of scattered photons

pixel_size_si = grid_spacing*binningsize;

photons_per_area = avg_photon_cnts_per_pix * epc / (qe * pixel_size_si^2);

number_photons = round(photons_per_area*absCross/alpha_fact);

%% 2D dist of scatterpoints 

dist_2D_conv_norm = recoil_blurring_dist(rg,grid_spacing,time,number_photons);

%% transfer func by fft

transfer_func = fft2(dist_2D_conv_norm);

end