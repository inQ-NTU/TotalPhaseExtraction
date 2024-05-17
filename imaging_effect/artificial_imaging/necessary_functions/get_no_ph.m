function nph = get_no_ph(grid_spacing_si,binningsize,avg_photon_cnts_per_pix,epc,qe,alpha_fact)

%% constants

absCross=2.906e-13; %m^2

%% calculate

pixel_size_si = grid_spacing_si*binningsize;
photons_per_area = avg_photon_cnts_per_pix * epc / (qe * pixel_size_si^2);
nph = photons_per_area*absCross/alpha_fact;

end