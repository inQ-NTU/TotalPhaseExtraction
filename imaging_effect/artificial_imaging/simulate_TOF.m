function [density2D, cloudwidths] = simulate_TOF(psi, TOF_si, grid_si, fringe_spacing, omega_trans, transversal_flag)

% Calculate 2D density of the gas after free expansion.
% psi has index structure [ realization number, condensate number , z-coordinate]
%
% transversal_flag.... 0 for VAndor images
%                      1 for TAndor
%
% density2D has index structure [z-coordinate, y-coordinate, realization number]

%% constants
hbar_si     = 6.626e-34/(2*pi);     
m_si        = 87*1.6605402e-27; % mass of Rb-87
as_si       = 98.98*52.917720859e-12; % s-wave scattering length


%% calculate transversal width of profile

a_perp      = sqrt(hbar_si/(m_si*omega_trans)); %harmonic oscillator ground state width
width_t     = a_perp*sqrt( 1 + (TOF_si*omega_trans)^2 ); % wavefunction width after free expansion


%% calculate density profiles

density2D   = zeros( length(grid_si), length(grid_si), size(psi,1) );
cloudwidths = zeros( length(grid_si), size(psi,1) );

for i = 1:size(psi, 1) % for each realization

    if size(psi,2) == 2 % double well
        % extract and propagate wavefunctions of each well individually
        psi_L       = remove_singleton( psi(i,1,:), 2 );
        psi_R       = remove_singleton( psi(i,2,:), 2 );

        psi_L_TOF   = func_free_TOF_propagation_1D( psi_L', grid_si, TOF_si );
        psi_R_TOF   = func_free_TOF_propagation_1D( psi_R', grid_si, TOF_si );

        %calculate 1D densities
        density_L   = conj(psi_L_TOF).*psi_L_TOF; % one dimensional density in left well
        density_R   = conj(psi_R_TOF).*psi_R_TOF; % --- in right well
        coh_TOF     = psi_L_TOF.*conj(psi_R_TOF); % interference term

        % correct transversal width of cloud for density-dependent expansion
        n_1D        = mean( [conj(psi_L).*psi_L ; conj(psi_R).*psi_R] )'; % in-trap density
        width       = width_t * ( 1 + 2*as_si*n_1D ).^(1/4);

        % half the distance between the wells
        halfd       = pi*m_si*a_perp^2 * width.^2 /(hbar_si*TOF_si*fringe_spacing); 

        % calculate 2D density profile
        if transversal_flag %for artificial TAndor images    
            dens1 = 1./(sqrt(pi)*width) .* exp(-(grid_si').^2 ./width.^2) * density_L; % 2D density profile of left well    
            dens2 = 1./(sqrt(pi)*width) .* exp(-(grid_si').^2 ./width.^2) * density_R; % ---- right well
            dens3 = zeros(size(dens2)); %interference term
        else %for artificial VAndor images
            dens1 = 1./(sqrt(pi)*width) .* exp(-(grid_si'+halfd).^2 ./width.^2) * density_L;
            dens2 = 1./(sqrt(pi)*width) .* exp(-(grid_si'-halfd).^2 ./width.^2) * density_R;
            dens3 = 2*real( 1./(sqrt(pi)*width) .* (exp(-(grid_si'.^2 + halfd.^2) ./width.^2) .* exp(1i*2*pi*grid_si'/fringe_spacing)) * coh_TOF); %interference term
        end

        density2D(:,:,i) = dens1 + dens2 + dens3;

    else % single well

        psi_L       = remove_singleton( psi(i,1,:), 2 );
        psi_TOF     = func_free_TOF_propagation_1D( psi_L', grid_si, TOF_si );
        density     = conj(psi_TOF).*psi_TOF; % one dimensional density 

        % correct transversal width of cloud for density-dependent expansion
        n_1D        = (conj(psi_L).*psi_L)'; % in-trap density
        width       = width_t * ( 1 + 2*as_si*n_1D ).^(1/4);

        density2D(:,:,i) = 1./(sqrt(pi)*width) .* exp(-(grid_si').^2 ./width.^2) * density; 
    end
    
    cloudwidths(:,i) = width;

end
    
end