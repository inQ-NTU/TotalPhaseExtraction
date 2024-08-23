%This is a code to stochastically sample phase and density fluctuation of
%1D Bose gases 
%The code is based on stochastic Bogoliubov sampling following the convention of
% reference below:
% Petrov, D. S. (2003). Bose-Einstein condensation in low-dimensional trapped gases (Doctoral dissertation, Universiteit van Amsterdam).

classdef class_bogoliubov_sampling < class_physical_parameters & handle
    properties
        temperature
        mean_density
        condensate_length
        coupling_strength_J
        wavevec_k
        g_coupling
        phase_samples
        density_fluct_samples
    end

    methods
        %1. Generate the constructor
        function obj = class_bogoliubov_sampling(temperature, mean_density_single_condensate, condensate_length, coupling_J)
            % Run the parent class constructor
            obj = obj@class_physical_parameters();
            % Set up some of the physics parameters, possibly to default
            if nargin < 2
                obj.mean_density = obj.max_longitudinal_density;
            else
                obj.mean_density = mean_density_single_condensate;
            end

            if nargin < 3
                obj.condensate_length = obj.default_condensate_length;
            else
                obj.condensate_length = condensate_length;
            end

            if nargin < 4
                obj.coupling_strength_J = obj.default_coupling_J;
            else
                obj.coupling_strength_J = coupling_J;
            end
            
            obj.temperature = temperature;
            obj.wavevec_k = 2*pi/obj.condensate_length;
            obj.g_coupling = 2*obj.hbar*obj.omega*obj.scattering_length;
        end


        %2. Stochastic sampling
        %Function to sample phase and density fluctuation - updated version
        function [phase_samples, density_fluct_samples] = generate_fluct_samples(obj, max_n_fourier, pixnumz, N_samples)
            if nargin < 4
                N_samples = 1;
            end
            %Define the grid and the sampling prefactors
            z_grid = linspace(-obj.condensate_length/2, obj.condensate_length/2, pixnumz);
            tunneling_energy = 2*obj.coupling_strength_J*obj.hbar;
            interaction_energy = 2*obj.g_coupling*obj.mean_density;

            %Initialize phase samples
            phase_samples = zeros(N_samples, pixnumz);
            density_fluct_samples = zeros(N_samples, pixnumz);
            %Start sampling
            for i = 1:N_samples
                for n = -max_n_fourier:max_n_fourier
                    if n~=0
                        kn = obj.wavevec_k*n;
                        Ekn = (obj.hbar*kn)^2/(2*obj.m);
                        bogoliubov_spectrum = sqrt((Ekn+tunneling_energy)*(Ekn+tunneling_energy+interaction_energy));
                        fkn_plus = sqrt(1/obj.condensate_length)*(bogoliubov_spectrum/(Ekn+tunneling_energy))^(1/2);
                        fkn_minus = sqrt(1/obj.condensate_length)*(bogoliubov_spectrum/(Ekn+tunneling_energy))^(-1/2);
                        thermal_occupation = (exp(bogoliubov_spectrum/(obj.kb*obj.temperature))-1)^(-1);
                        
                        z0 = randn; z1 = randn; 
                        bn = sqrt(thermal_occupation/2)*(z0+1j*z1);
                        
                        for j = 1:pixnumz
                            phase_term = sqrt(1/(4*obj.mean_density))*fkn_plus*bn*exp(1j*kn*z_grid(j));
                            density_term = 1j*sqrt(obj.mean_density)*fkn_minus*bn*exp(1j*kn*z_grid(j));
                            phase_samples(i,j) = phase_samples(i,j)+phase_term+conj(phase_term);
                            density_fluct_samples(i,j) = density_fluct_samples(i,j)+density_term+conj(density_term);
                        end
                    end
                end
            end
            obj.phase_samples = phase_samples;
            obj.density_fluct_samples = density_fluct_samples;
        end

        %3. A function to coarse grain/interpolate a given phase sample
        function coarse_fine_profile = coarse_fine_grain(obj, target_longitudinal_points, profile)
            if nargin < 3
                profile = obj.phase_samples;
            end
            nmb_longitudinal_points = size(profile,2);
            [Z] = ndgrid(1:nmb_longitudinal_points)';
            nmb_samples = size(profile,1);
            coarse_fine_profile = zeros(nmb_samples,target_longitudinal_points);
            for i = 1:nmb_samples
                profile_interpolant = griddedInterpolant(Z,profile(i,:));
                coarse_fine_profile(i,:) = profile_interpolant(linspace(1,nmb_longitudinal_points, target_longitudinal_points));
            end
        end

        %4. A function to analytically compute correlation function
        function corr = corr_function(obj, max_fourier, pixnumz)
            z_grid = linspace(-obj.condensate_length/2, obj.condensate_length/2, pixnumz);
            interaction_energy = 2*obj.g_coupling*obj.mean_density;
            phase_var = zeros(1, pixnumz);
            for n = -max_fourier:max_fourier
                if n~= 0
                    kn = obj.wavevec_k*n;
                    Ekn = (obj.hbar*kn)^2/(2*obj.m);
                    bogoliubov_spectrum = sqrt(Ekn*(Ekn+interaction_energy));
                    thermal_occupation = (exp(bogoliubov_spectrum/(obj.kb*obj.temperature))-1)^(-1);
                    phase_var = phase_var + (bogoliubov_spectrum/Ekn)*(1-cos(kn.*z_grid))*(1+2*thermal_occupation);
                end
            end
            prefactor = (2*obj.mean_density*obj.condensate_length)^(-1);
            phase_var = phase_var*prefactor;
            corr = exp(-phase_var/2);
        end
    end %end of methods
end %end of class