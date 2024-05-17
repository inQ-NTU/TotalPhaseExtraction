function atoms = func_recoil_blurring(atoms,time,avg_photon_cnts_per_pix,epc,qe,rg,grid_spacing,binningsize,alpha_fact)

transfer_func = calc_recoil_blurring_transfer_func(time,avg_photon_cnts_per_pix,epc,qe,rg,grid_spacing,binningsize,alpha_fact);

atoms = ifftshift(ifft2(transfer_func.*fft2(atoms)));

end