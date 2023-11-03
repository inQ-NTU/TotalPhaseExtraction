%This is a class to stochastically sample phase profiles in 1D Bose gas
%interferometry experiment

%The scheme is based on the reference below
%Stimming, H. P., Mauser, N. J., Schmiedmayer, J., & Mazets, I. E. (2010). 
% Fluctuations and stochastic processes in one-dimensional many-body quantum systems. PRL , 105(1), 015301.

classdef class_bogoliubov_phase_sampling < class_physical_parameters & handle
    properties
        temperature
        mean_density
        condensate_length
        coupling_strength_J

        fourier_cosine_coeffs
        fourier_sine_coeffs
        phase_samples
    end

    methods
        %Generate the constructor
        function obj = class_bogoliubov_phase_sampling(temperature, coupling_J, mean_density, condensate_length)
            % Run the parent class constructor
            obj = obj@class_physical_parameters();
            % Set up some of the physics parameters, possibly to default
            if nargin < 2
                obj.coupling_strength_J = obj.default_coupling_J;
            else
                obj.coupling_strength_J = coupling_J;
            end
            
            if nargin < 3
                obj.mean_density = obj.max_longitudinal_density;
            else
                obj.mean_density = mean_density;
            end

            if nargin < 4
                obj.condensate_length = obj.default_condensate_length;
            else
                obj.condensate_length = condensate_length;
            end
            obj.temperature = temperature;
        end


        %Stochastic sampling 
        function phase_samples = generate_samples(obj, max_n_fourier, z_resolution, N_samples)
            if nargin < 4
                N_samples = 1;
            end
            z_grid = linspace(-obj.condensate_length/2, obj.condensate_length/2, z_resolution);
            h = obj.hbar*2*pi; 
            prefactor = sqrt(4*obj.m*obj.kb*obj.temperature*obj.condensate_length/(obj.mean_density*h^2));
            tunnel_coupling_term = 2*obj.m*(obj.condensate_length)^2*obj.coupling_strength_J/(pi*h);
            phase_samples = zeros(N_samples,z_resolution);
            cosine_coeffs = zeros(N_samples, max_n_fourier);
            sine_coeffs = zeros(N_samples, max_n_fourier);
            for j = 1:N_samples
                for n = 1:max_n_fourier
                    random_num1 = abs(log(rand()));
                    random_num2 = rand();
                    denom = n^2 + tunnel_coupling_term;
                    cosine_coeffs(j,n) = prefactor*sqrt(random_num1/denom)*sin(2*pi*random_num2);
                    sine_coeffs(j,n) = prefactor*sqrt(random_num1/denom)*cos(2*pi*random_num2);
                    for i = 1:z_resolution
                        phase_samples(j,i) = phase_samples(j,i) + sqrt(random_num1/denom)*sin(2*pi*n*z_grid(i)/obj.condensate_length+2*pi*random_num2);
                    end
                end
            end
            phase_samples = phase_samples*prefactor; %Multiplying temperature dependent prefactor
            obj.phase_samples = phase_samples;
            obj.fourier_cosine_coeffs = cosine_coeffs;
            obj.fourier_sine_coeffs = sine_coeffs;
        end
    end %end of methods
end %end of class