function ctf2D = ctf_gauss_cloud(nutg,Rz,dz,numap)

% Rz.....Gauss sigma of the wavefunction, not of the density!!!!!!!!!

lambda = 780e-9;

ctf2D = exp(-(pi*lambda/2*nutg.^2*Rz).^2).*exp(-1i*pi*lambda*nutg.^2*dz).*heaviside(numap/lambda - nutg);

end