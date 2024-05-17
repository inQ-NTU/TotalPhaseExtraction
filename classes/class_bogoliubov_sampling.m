%This is a class to stochastically sample phase profiles in 1D Bose gas
%interferometry experiment

%The scheme is based on the reference below
%Stimming, H. P., Mauser, N. J., Schmiedmayer, J., & Mazets, I. E. (2010). 
% Fluctuations and stochastic processes in one-dimensional many-body quantum systems. PRL , 105(1), 015301.

classdef class_bogoliubov_sampling < class_physical_parameters & handle
    properties
        temperature
        mean_density
        condensate_length
        coupling_strength_J
        wavevec_k

        %fourier_cosine_coeffs
        %fourier_sine_coeffs
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
        end


        %2. Stochastic sampling

        %Function to sample phase - updated version
        function phase_samples = generate_phase_samples(obj, max_n_fourier, pixnumz, N_samples)
            if nargin < 4
                N_samples = 1;
            end
            
            %Define the grid and the sampling prefactors
            z_grid = linspace(-obj.condensate_length/2, obj.condensate_length/2, pixnumz);
            overall_prefactor = sqrt(2*obj.kb*obj.temperature/(obj.mean_density*obj.condensate_length));
            tunneling_energy = 2*obj.coupling_strength_J*obj.hbar;
            %interaction_energy = 2*obj.g_coupling*obj.mean_density;

            %Initialize phase samples
            phase_samples = zeros(N_samples, pixnumz);
            
            %Start sampling
            for i = 1:N_samples
                for n = -max_n_fourier:max_n_fourier
                    if n~=0
                        kn = obj.wavevec_k*n;
                        Ekn = (obj.hbar*kn)^2/(2*obj.m)+tunneling_energy;
                        random_num1 = abs(log(rand()));
                        random_num2 = rand();
                        for j = 1:pixnumz
                            phase_samples(i,j) = phase_samples(i,j)+sqrt(random_num1/Ekn)*sin(kn*z_grid(j)+2*pi*random_num2);
                        end
                    end
                end
            end
            phase_samples = overall_prefactor*phase_samples;
            obj.phase_samples = phase_samples;
        end

        %Function to sample density fluctuation
        function density_fluct_samples = generate_density_fluct_samples(obj, max_n_fourier, pixnumz, N_samples)
            if nargin < 4
                N_samples = 1;
            end

            %Define the grid and the sampling prefactors
            z_grid = linspace(-obj.condensate_length/2, obj.condensate_length/2, pixnumz);
            overall_prefactor = sqrt(2*obj.kb*obj.temperature/(obj.mean_density*obj.condensate_length));
            tunneling_energy = 2*obj.coupling_strength_J*obj.hbar;
            g_coupling = 2*obj.hbar*obj.omega*obj.scattering_length;
            interaction_energy = 2*g_coupling*obj.mean_density;

            %Initialize density fluctuation samples
            density_fluct_samples = zeros(N_samples, pixnumz);
            
            %Start sampling
            for i = 1:N_samples
                for n = -max_n_fourier:max_n_fourier
                    if n~=0
                        kn = obj.wavevec_k*n;
                        Ekn = ((obj.hbar*kn)^2)/(2*obj.m)+tunneling_energy+interaction_energy;
                        random_num1 = abs(log(rand()));
                        random_num2 = rand();
                        for j = 1:pixnumz
                            density_fluct_samples(i,j) = density_fluct_samples(i,j)+sqrt(random_num1/Ekn)*sin(kn*z_grid(j)+2*pi*random_num2);
                        end
                    end
                end
            end
            density_fluct_samples = overall_prefactor*obj.mean_density*density_fluct_samples;
            obj.phase_samples = density_fluct_samples;
        end
        %Version 15 May 2024
        %function phase_samples = generate_phase_samples(obj, max_n_fourier, z_resolution, N_samples)
        %    if nargin < 4
        %        N_samples = 1;
        %    end
        %    z_grid = linspace(-obj.condensate_length/2, obj.condensate_length/2, z_resolution);
        %    h = obj.hbar*2*pi; 
        %    prefactor = sqrt(4*obj.m*obj.kb*obj.temperature*obj.condensate_length/(obj.mean_density*h^2));
        %    tunnel_coupling_term = 2*obj.m*(obj.condensate_length)^2*obj.coupling_strength_J/(pi*h);
        %    phase_samples = zeros(N_samples,z_resolution);
        %    cosine_coeffs = zeros(N_samples, 2*max_n_fourier);
        %   sine_coeffs = zeros(N_samples, 2*max_n_fourier);
        %    for j = 1:N_samples
        %        for n = -max_n_fourier:max_n_fourier
        %            if n ~= 0
        %                random_num1 = abs(log(rand()));
        %                random_num2 = rand();
        %                denom = n^2 + tunnel_coupling_term;
        %                if n<0
        %                    cosine_coeffs(j,n+max_n_fourier+1) = prefactor*sqrt(random_num1/denom)*sin(2*pi*random_num2);
        %                    sine_coeffs(j,n+max_n_fourier+1) = prefactor*sqrt(random_num1/denom)*cos(2*pi*random_num2);
        %                elseif n>0
        %                    cosine_coeffs(j,n+max_n_fourier) = prefactor*sqrt(random_num1/denom)*sin(2*pi*random_num2);
        %                    sine_coeffs(j,n+max_n_fourier) = prefactor*sqrt(random_num1/denom)*cos(2*pi*random_num2);
        %                end
        %                for i = 1:z_resolution
        %                    phase_samples(j,i) = phase_samples(j,i) + sqrt(random_num1/denom)*sin(2*pi*n*z_grid(i)/obj.condensate_length+2*pi*random_num2);
        %                end
        %            end
        %        end
        %    end
        %    phase_samples = phase_samples*prefactor; %Multiplying temperature dependent prefactor
        %    obj.phase_samples = phase_samples;
        %    obj.fourier_cosine_coeffs = cosine_coeffs;
        %    obj.fourier_sine_coeffs = sine_coeffs;
        %end

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