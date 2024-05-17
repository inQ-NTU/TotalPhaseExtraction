function [tf2D,sigma_rec] = get_tf_recoil(input_data)

% user input
grid_si = input_data.grid_si;
binningsize = input_data.binningsize;
avg_photon_cnts_per_pix = input_data.avg_photon_cnts_per_pix;
epc = input_data.epc;
qe = input_data.qe;
alpha_fact = input_data.alpha_fact;
n_sub = input_data.no_push_subdivision;    
t_imag = input_data.imaging_time;

%% calculate some helpful quantities

no_gridpoints = length(grid_si);
grid_spacing_si = grid_si(2) - grid_si(1);
extend = no_gridpoints * grid_spacing_si;

%% calculate number of scattered photons per atom

nph = get_no_ph(grid_spacing_si,binningsize,avg_photon_cnts_per_pix,epc,qe,alpha_fact);

%% claculate the recoil transfer function

[tf2D,sigma_rec] = tf_recoil(t_imag,nph,n_sub,extend,no_gridpoints);

end