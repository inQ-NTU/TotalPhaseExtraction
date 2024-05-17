function [psi, phi, drho] = generate_equilibrium_wavefunctions( N_samples, ...      % number of realisations
                                                                dens_profile, ...   % background density profile
                                                                grid_si, ...        % spatial grid
                                                                T_si, ...           % temperature
                                                                J_si, ...           % tunnel couplings in Hz
                                                                periodic_BC, ...    % flag for periodic boundary conditions 
                                                                omegaT_si, ...      % transverse trapping frequency
                                                                N_wells )           % number of wells (typically 1 or 2)
                                                            
% Given a background density profile, generate N_samples realisations
% of wavefunctions with density and phase fluctuations given by the
% equilibrium theory.
    

%% 

N_gridpoints    = length(grid_si);
g1d             = calc_g0(omegaT_si);
grid_spacing    = grid_si(2) - grid_si(1);

% make sure density profile is row vectors
dens_profile    = reshape(dens_profile, 1, []);


%% Generate equilibrium wave function (psi)

% ignore regions with small densities to calculate fluctuations
% (theory not valid in those regimes)
idx_bogo    = dens_profile > 0.1*max(dens_profile); 
profile_bogo= dens_profile(idx_bogo);

% sample density fluctuations using Bogoliubov theory    
dens_fluct_init_si  = bogo_density_fluct_shots( T_si, ...
                                                profile_bogo, ...
                                                g1d*ones(size(profile_bogo)), ...
                                                J_si, ...
                                                grid_spacing, ...
                                                periodic_BC, ...
                                                N_samples*N_wells);


dens_fluct = zeros(N_samples*N_wells, N_gridpoints);
dens_fluct(:,idx_bogo) = dens_fluct_init_si; 


% sample phase fluctuations and construct wavefunction psi
psi = 1i*ones(N_samples*N_wells, N_gridpoints);
phi = zeros(N_samples*N_wells, N_gridpoints);

for j = 1:N_samples*N_wells
    % Ornstein-Uhlenbeck simulation of phase distribution
    phase_distribution = DensityDependantOU_SingleQuasiBECFun(  T_si, ...
                                                                N_gridpoints, ...
                                                                grid_spacing, ...
                                                                dens_profile); 

    profile_temp = abs(dens_profile + squeeze(dens_fluct(j,:)));

    % psi = sqrt(\rho + \delta\rho) exp(i\phi)
    psi(j, :)   = sqrt(profile_temp) .* exp(1i*phase_distribution');
    phi(j,:)    = phase_distribution;
end

% reshape psi to correct index structure
psi     = reshape(psi, [N_samples, N_wells, N_gridpoints]);
phi     = reshape(phi, [N_samples, N_wells, N_gridpoints]);
drho    = reshape(dens_fluct, [N_samples, N_wells, N_gridpoints]);


end 
