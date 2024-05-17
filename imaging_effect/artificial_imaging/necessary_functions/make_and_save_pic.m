function [finalP_atoms,finalP_noatoms] = make_and_save_pic(psi_L,psi_R,TOF_si,grid_si,fringeSpacing,omega_trans,transversal_flag,...
    recoil_flag,avg_photon_cnts_per_pix,epc,qe,binningsize,alpha_fact,shotnoise_flag,ctf2D,tf_recoil,imaging_shift,shift_dir,save_path,real_ind,save_pics_flag,background_cnts)

no_gridpoints = length(grid_si);

[atoms,~] = func_create_atomic_density(psi_L,psi_R,TOF_si,grid_si,fringeSpacing,omega_trans,transversal_flag);

if ~recoil_flag
    [finalP_atoms,finalP_noatoms] = func_create_art_image_non_sat...
        (atoms,no_gridpoints,avg_photon_cnts_per_pix,shotnoise_flag,binningsize,ctf2D,epc,qe,alpha_fact,imaging_shift,shift_dir);
else
    [finalP_atoms,finalP_noatoms] = func_create_art_image_non_sat...
        (atoms,no_gridpoints,avg_photon_cnts_per_pix,shotnoise_flag,binningsize,ctf2D,epc,qe,alpha_fact,imaging_shift,shift_dir,tf_recoil);
end

% % just for testing
% absCross=2.906e-13; %m^2
% absPic_show=-(log(finalP_atoms./(finalP_noatoms+eps)));
% atomnum_test = alpha_fact*(2e-6)^2*sum(sum(absPic_show))/absCross;
% atomnum_test_density = (0.5e-6)^2*sum(sum(atoms));
% atomnum_test_should_be = (sum(abs(psi_L(real_ind,:)).^2) + sum(abs(psi_R(real_ind,:)).^2))*0.5e-6;

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