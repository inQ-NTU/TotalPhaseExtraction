function [atom_density_2D,finalP_atoms,finalP_noatoms] = make_and_analyze_pic(psi_L,psi_R,TOF_si,grid_si,fringeSpacing,omega_trans,transversal_flag,...
    recoil_flag,avg_photon_cnts_per_pix,epc,qe,binningsize,alpha_fact,shotnoise_flag,ctf2D,tf_recoil,imaging_shift,shift_dir)

%% constants

absCross=2.906e-13; %m^2

%%

grid_spacing_si = grid_si(2) - grid_si(1);
no_gridpoints = length(grid_si);

[atoms,~] = func_create_atomic_density(psi_L,psi_R,TOF_si,grid_si,fringeSpacing,omega_trans,transversal_flag);

if ~recoil_flag
    [finalP_atoms,finalP_noatoms] = func_create_art_image_non_sat...
        (atoms,no_gridpoints,avg_photon_cnts_per_pix,shotnoise_flag,binningsize,ctf2D,epc,qe,alpha_fact,imaging_shift,shift_dir);
else
    [finalP_atoms,finalP_noatoms] = func_create_art_image_non_sat...
        (atoms,no_gridpoints,avg_photon_cnts_per_pix,shotnoise_flag,binningsize,ctf2D,epc,qe,alpha_fact,imaging_shift,shift_dir,tf_recoil);
end

absPic_show=-(log(finalP_atoms./(finalP_noatoms+eps)));
atom_density_2D = alpha_fact*(grid_spacing_si*binningsize)^2*absPic_show/absCross;

end