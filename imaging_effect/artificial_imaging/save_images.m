function save_images(save_path, atompics, backpics)

% atompics and backpics are 3D arrays with dimensions
% [N_xpixels, N_ypixels, N_pictures]

for i = 1:size(atompics)
    % save artificial images
    pic_counts_atomcloud    = uint16(atompics(:,:,i));
    name_atomcloud          = [save_path, num2str(real_ind), '-atomcloud.tif'];
    imwrite_no_overwrite( pic_counts_atomcloud, name_atomcloud );

    pic_counts_withoutatoms = uint16(backpics(:,:,i));
    name_withoutatoms       = [save_path, num2str(real_ind), '-withoutatoms.tif'];
    imwrite_no_overwrite( pic_counts_withoutatoms, name_withoutatoms );
end