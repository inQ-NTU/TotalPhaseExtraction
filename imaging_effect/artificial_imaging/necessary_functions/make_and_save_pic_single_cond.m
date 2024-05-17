function [finalP_atoms,finalP_noatoms] = make_and_save_pic_single_cond(psi,TOF_si,grid_si,omega_trans,...
    recoil_flag,avg_photon_cnts_per_pix,epc,qe,binningsize,alpha_fact,shotnoise_flag,ctf2D,tf_recoil,imaging_shift,shift_dir,save_path,real_ind,save_pics_flag,background_cnts)

no_gridpoints = length(grid_si);

[atoms,~] = func_create_atomic_density_single_cond(psi,TOF_si,grid_si,omega_trans);

%---- just for testing
% test1 = sum(atoms,1)*0.5e-6;
% test2 = abs(psi).^2;

% figure(1)
% clf
% hold on
% plot(test1)
% plot(test2,'--')


if ~recoil_flag
    [finalP_atoms,finalP_noatoms] = func_create_art_image_non_sat...
        (atoms,no_gridpoints,avg_photon_cnts_per_pix,shotnoise_flag,binningsize,ctf2D,epc,qe,alpha_fact,imaging_shift,shift_dir);
else
    [finalP_atoms,finalP_noatoms] = func_create_art_image_non_sat...
        (atoms,no_gridpoints,avg_photon_cnts_per_pix,shotnoise_flag,binningsize,ctf2D,epc,qe,alpha_fact,imaging_shift,shift_dir,tf_recoil);
end

if save_pics_flag
    %save artificial images, the background counts are added to each pixel
    
    pic_counts_atomcloud = uint16(finalP_atoms')+background_cnts;
    name_atomcloud =[save_path,num2str(real_ind),'-atomcloud.tif'];
    imwrite_no_overwrite(pic_counts_atomcloud,name_atomcloud);
    
    pic_counts_withoutatoms = uint16(finalP_noatoms')+background_cnts;
    name_withoutatoms =[save_path,num2str(real_ind),'-withoutatoms.tif'];
    imwrite_no_overwrite(pic_counts_withoutatoms,name_withoutatoms);
end

end