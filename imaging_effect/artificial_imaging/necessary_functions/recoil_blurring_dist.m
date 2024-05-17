function dist_2D_conv_norm = recoil_blurring_dist(rg,grid_spacing,time,number_photons)

%% constants

vrec = 5.8845e-3;

%% 2D dist of scatterpoints 

%to avoid inf
rg(rg == 0) = grid_spacing;

c = 1/9 * vrec^2 * (time/number_photons)^2;
c2 = rg.^2/(2*c);
c1 = 1/number_photons*(2*pi*c)^(-1);
dist_2D = c1*gamma(2/3).*gammainc(c2/number_photons^3,2/3,'upper')./(3*c2.^(2/3));

dist_2D_conv_norm = dist_2D/sum(sum(dist_2D)); 
% note: normalizing by multiplying with gridspacing^2 usually does not work 
%as the function is quite narrow compared to the grid_spacing

end