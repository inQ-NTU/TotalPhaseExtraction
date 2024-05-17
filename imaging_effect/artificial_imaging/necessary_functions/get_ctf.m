function [ctf2D,rg] = get_ctf(input_data)

% user input that is always needed
grid_si = input_data.grid_si;
use_gaussian_psf = input_data.use_gaussian_psf;

no_gridpoints = length(grid_si);
grid_spacing_si = grid_si(2) - grid_si(1);

%% calculate

% make a grid
[xg,yg]=meshgrid(grid_si,grid_si);
rg = sqrt(xg.^2 + yg.^2);

if use_gaussian_psf
    % user input that is only needed if using a gaussian psf
    psf_sigma = input_data.psf_sigma;
    
    psf2D = exp( -(rg.^2)./(2*psf_sigma^2));
    psf2D = psf2D/sum(sum(psf2D));
    ctf2D = fftshift(fft2(fftshift(psf2D))); %no idea why one needs fftshift(fft2(fftshift(...))), figured out by trail and error
else   
    % user input that is only needed if not using a gaussian psf
    cloud_sigma_for_psf = input_data.cloud_sigma_for_psf;
    TOF_si = input_data.TOF_si;
    omega_trans = input_data.omega_trans;
    defocus_um = input_data.defocus_um;
    numap = input_data.numap;
    push_flag = input_data.push_flag;
    
    %--- create grid in Fourier-Space
    size_of_pic = no_gridpoints*grid_spacing_si;
    nu_spacing = 1/size_of_pic;
    if mod(no_gridpoints,2)
        nu_grid = get_symmetric_z_axis(no_gridpoints,nu_spacing);
    else
        nu_grid = -(no_gridpoints/2)*nu_spacing:nu_spacing:(no_gridpoints/2-1)*nu_spacing;
    end
    [nuxg,nuyg]=meshgrid(nu_grid,nu_grid);
    nutg = sqrt(nuxg.^2 + nuyg.^2);
    
    if isnan(cloud_sigma_for_psf)
        cloud_sigma_for_psf = func_ideal_TOF_width(TOF_si,omega_trans);
    end
    
    if ~push_flag
        ctf2D = ctf_gauss_cloud(nutg,cloud_sigma_for_psf,defocus_um*1e-6,numap);
    else
        %user input that is only needed when considering the imaging push
        imaging_time = input_data.imaging_time;
        no_push_subdivision = input_data.no_push_subdivision;
        avg_photon_cnts_per_pix = input_data.avg_photon_cnts_per_pix;
        epc = input_data.epc;
        qe = input_data.qe;
        alpha_fact = input_data.alpha_fact;
        binningsize = input_data.binningsize;
        
        %for considering that the cloud is pushed through focus
        defocus_arr_si = imaging_push(imaging_time,no_push_subdivision,avg_photon_cnts_per_pix,epc,qe,alpha_fact,grid_spacing_si,binningsize);
        defocus_arr_si = defocus_arr_si + defocus_um*1e-6;
        
        ctf2D = zeros(no_gridpoints,no_gridpoints,no_push_subdivision);
        for push_ind=1:no_push_subdivision
            ctf2D(:,:,push_ind) = ctf_gauss_cloud(nutg,cloud_sigma_for_psf,defocus_arr_si(push_ind),numap);
        end
    end
end

end