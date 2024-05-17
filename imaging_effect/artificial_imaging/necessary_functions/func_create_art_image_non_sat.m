function [finalP_atoms,finalP_noatoms] = func_create_art_image_non_sat...
    (atoms,no_gridpoints,avg_cnts_per_pix,shotnoise_flag,binningsize,CTF_func2D,epc,qe,alpha_fact,imaging_shift,shift_dir,ctf_recoil_blurr)
%%
% 
% input:
%   imaging_shift.....shift in gridpoints of the atoms, not in pixels
%   shift_dir.....can be 1 or 2, for Tandor 2 is the proper value

if nargin < 10
    imaging_shift = 0;
end

if nargin < 12
    recoil_flag = 0;
else
    recoil_flag = 1;
end

%% constants

absCross=2.906e-13; %m^2

%%  create artificial image

avgphotons = avg_cnts_per_pix * epc / (qe * binningsize^2);

%create the beam profile
beam=ones(no_gridpoints,no_gridpoints);

% to have the possibility to consider that the cloud is pushed out of focus
no_push_subdivision = size(CTF_func2D,3);

finalP_atomsPSFI_arr = zeros(no_push_subdivision,no_gridpoints,no_gridpoints);

if imaging_shift > 1
    shift_arr = linspace(0,imaging_shift,no_push_subdivision);
    shift_arr = shift_arr -  mean(shift_arr);
    shift_arr = round(shift_arr);
end

for push_ind = 1:no_push_subdivision
    CTF_func2D_to_use = squeeze(CTF_func2D(:,:,push_ind));
    
    if recoil_flag
        ctf_recoil_to_use = squeeze(ctf_recoil_blurr(push_ind,:,:));
        atoms_temp = ifft2(ctf_recoil_to_use.*fft2(atoms));
    else
        atoms_temp = atoms;
    end
    
    %put atoms on beam
    finalP_atoms = beam.*exp(-absCross/alpha_fact * atoms_temp/2); % factor 1/2 because it's the field, not the intensity
    
    %apply PSF/CTF
    finalP_atomsPSF = ifft2(CTF_func2D_to_use.*fft2(finalP_atoms));
    
    %calculate intensity
    finalP_atomsPSFI_arr(push_ind,:,:) = abs(finalP_atomsPSF).^2 * avgphotons;
    
    if imaging_shift > 0
        finalP_atomsPSFI_arr(push_ind,:,:) = circshift(finalP_atomsPSFI_arr(push_ind,:,:),shift_arr(push_ind),shift_dir+1);
    end
end

finalP_atomsPSFI = squeeze(mean(finalP_atomsPSFI_arr,1));

finalP_noatomsPSFI = beam * avgphotons; %background image = beam

%bin the image
finalP_atomsPSFI_binned = binMatrix(finalP_atomsPSFI,binningsize);
finalP_noatomsPSFI_binned = binMatrix(finalP_noatomsPSFI,binningsize);

%consider quantum efficiency
finalP_atomsPSFI_binned = finalP_atomsPSFI_binned*qe;
finalP_noatomsPSFI_binned = finalP_noatomsPSFI_binned*qe;

%add shotnoise
finalP_atomsPSFI_binned_noise = double(poissrnd(finalP_atomsPSFI_binned));
finalP_noatomsPSFI_binned_noise = double(poissrnd(finalP_noatomsPSFI_binned));

%choose whether to consider noise or no
if(shotnoise_flag)
    finalP_atoms = finalP_atomsPSFI_binned_noise/epc;
    finalP_noatoms = finalP_noatomsPSFI_binned_noise/epc;
else
    finalP_atoms = finalP_atomsPSFI_binned/epc;
    finalP_noatoms = finalP_noatomsPSFI_binned/epc;
end

end