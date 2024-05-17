function [tf2D,sigma_rec] = tf_recoil(t_imag,nph,n_sub,x_extend,no_gridpoints)
%% constants

vrec = 5.8845e-3; %recoil velocity

%% division into sub-pics

timestep = t_imag/n_sub; 
t_sub = timestep/2:timestep:t_imag;

nph_step = nph/n_sub;
nph_sub = nph_step/2:nph_step:nph;

%% calculate ctf

% make a grid in k-space
k_spacing = 2*pi/x_extend;
k_grid_1D = get_fft_grid(k_spacing,no_gridpoints);
[kxg,kyg]=meshgrid(k_grid_1D,k_grid_1D);
krg = sqrt(kxg.^2 + kyg.^2);

sigma_rec = zeros(1,n_sub);
tf2D = zeros(n_sub,no_gridpoints,no_gridpoints);
for i = 1:n_sub
    sigma_rec(i) = sqrt(nph_sub(i)*vrec^2*t_sub(i)^2/9);
    ctf2D_temp = exp(-krg.^2 * sigma_rec(i)^2 / 2);
    ctf2D_temp = ifftshift(ctf2D_temp);
    tf2D(i,:,:) = ctf2D_temp;
end

end