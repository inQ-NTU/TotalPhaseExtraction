%This is a class to stochastically sample phase profiles in 1D Bose gas
%interferometry experiment

%The scheme is based on the references below:
%Stimming, H. P., Mauser, N. J., Schmiedmayer, J., & Mazets, I. E. (2010). 
%Fluctuations and stochastic processes in one-dimensional many-body quantum systems. PRL , 105(1), 015301.
%Stimming, H. P., Mauser, N. J., Schmiedmayer, J., & Mazets, I. E. (2011). 
% Dephasing in coherently split quasicondensates. PRA, 83(2), 023618.

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
        function obj = class_bogoliubov_sampling(temperature, mean_density, condensate_length, coupling_J)
            % Run the parent class constructor
            obj = obj@class_physical_parameters();
            % Set up some of the physics parameters, possibly to default
            if nargin < 2
                obj.mean_density = obj.max_longitudinal_density;
            else
                obj.mean_density = mean_density;
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
            overall_prefactor_phase = sqrt(2*obj.kb*obj.temperature/(obj.mean_density*obj.condensate_length));
            overall_prefactor_density = sqrt(2*obj.kb*obj.temperature*obj.mean_density/obj.condensate_length);
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
                        Ekn_phase = (obj.hbar*kn)^2/(2*obj.m)+tunneling_energy;
                        Ekn_density = (obj.hbar*kn)^2/(2*obj.m)+tunneling_energy+interaction_energy;
                        random_num1 = abs(log(rand()));
                        random_num2 = rand();
                        for j = 1:pixnumz
                            phase_samples(i,j) = phase_samples(i,j)+sqrt(random_num1/Ekn_phase)*sin(kn*z_grid(j)+2*pi*random_num2);
                            density_fluct_samples(i,j) = density_fluct_samples(i,j)+sqrt(random_num1/Ekn_density)*cos(kn*z_grid(j)+2*pi*random_num2);
                        end
                    end
                end
            end
            phase_samples = overall_prefactor_phase*phase_samples;
            density_fluct_samples = overall_prefactor_density*density_fluct_samples;
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
    end %end of methods
end %end of class