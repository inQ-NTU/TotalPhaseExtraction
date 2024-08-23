classdef class_common_phase_extraction < class_physical_parameters & handle
    properties
        %input
        density_ripple
        mean_density
        expansion_time
        z_grid
        boundary_cut_ratio
        boundary_condition
        well_number
        
        %inferred intermediate variables
        mean_density_gradient
        condensate_length
        wavevec_k
        lt %length scale of longitudinal expansion
        dimensionless_grid
        dimensionless_mean_density_gradient
        dimensionless_density_ripple

        %output
        common_phase_profile
    end

    methods
        %1. Declare the constructor
        function obj = class_common_phase_extraction(density_ripple, mean_density, z_grid, expansion_time, boundary_condition, well_nmb)
            if nargin < 4
                obj.expansion_time = obj.default_expansion_time;
            else
                obj.expansion_time = expansion_time;
            end

            if nargin < 5
                obj.boundary_condition = 'Periodic';
            else
                obj.boundary_condition = boundary_condition;
            end

            if nargin < 6
                obj.well_number = 2;
                obj.lt = sqrt(obj.hbar*obj.expansion_time/(2*obj.m)); 
            else
                obj.well_number = well_nmb;
                obj.lt = sqrt(obj.hbar*obj.expansion_time/(well_nmb*obj.m));
            end

            obj.mean_density = mean_density;
            %obj.z_grid = z_grid;
            obj.condensate_length = abs(z_grid(end)-z_grid(1));
            obj.mean_density_gradient = gradient(obj.mean_density);

            %Setting up boundary condition
            if strcmp(obj.boundary_condition,'Periodic')
                obj.wavevec_k = 2*pi/obj.condensate_length;
                if z_grid(end) < obj.condensate_length/2
                    shift = obj.condensate_length/2-z_grid(end);
                    obj.z_grid = z_grid + shift;
                elseif z_grid(end) > obj.condensate_length
                    shift = z_grid(end) - obj.condensate_length/2;
                    obj.z_grid = z_grid - shift;
                else
                    obj.z_grid = z_grid;
                end
            elseif strcmp(obj.boundary_condition, 'Neumann')| strcmp(obj.boundary_condition, 'Closed')
                obj.wavevec_k = pi/obj.condensate_length;
                obj.z_grid = z_grid - z_grid(1);
            else 
                disp(['Sorry, input boundary condition is unknown, ...' ...
                    'using the default periodic boundary instead'])
                obj.wavevec_k = 2*pi/obj.condensate_length;
                if z_grid(end) < obj.condensate_length/2
                    shift = obj.condensate_length/2-z_grid(end);
                    obj.z_grid = z_grid + shift;
                elseif z_grid(end) > obj.condensate_length
                    shift = z_grid(end) - obj.condensate_length/2;
                    obj.z_grid = z_grid - shift;
                end
            end

            obj.density_ripple = density_ripple;            
            obj.dimensionless_grid = obj.z_grid/obj.lt;
            obj.dimensionless_mean_density_gradient = gradient(obj.mean_density, obj.dimensionless_grid)./obj.mean_density;
            
            dimensionless_density_ripple = 1-(obj.density_ripple./obj.mean_density);
            dz = obj.z_grid(2) - obj.z_grid(1);
            obj.dimensionless_density_ripple = dimensionless_density_ripple - (dz/obj.condensate_length)*sum(dimensionless_density_ripple);
            
        end

        %2. Function to extract common phase spectrum when we ignore
        %gradient of mean density
        function [cosineCoeffs, sineCoeffs] = extract_com_spectrum(obj, ext_cutoff)

            %initialize the outputs
            cosineCoeffs = zeros(1, ext_cutoff);
            sineCoeffs = zeros(1, ext_cutoff);
            complex_fourier_coeffs = zeros(1,ext_cutoff);

            %Compute and initialize parameters
            dz = obj.z_grid(2) - obj.z_grid(1);
            ripple = obj.dimensionless_density_ripple;
            
            %Making sure zero constant shift in Fourier decomposition
            ripple = ripple - sum(ripple)*dz/obj.condensate_length; 
            
            %Start computing Fourier coefficients
            for n = 1:ext_cutoff
                fn = 0;
                for i = 1:length(obj.z_grid)
                    fn = fn + ripple(i)*exp(1j*n*obj.wavevec_k*obj.z_grid(i));
                end
                complex_fourier_coeffs(n) = (dz/obj.condensate_length)*fn;
            end

            %Assigning outputs
            for n = 1:ext_cutoff
                complex_fourier_coeffs(n) = (-2*complex_fourier_coeffs(n))/((obj.wavevec_k*n*obj.lt)^2);
                if strcmp(obj.boundary_condition, 'Neumann')
                    cosineCoeffs(n) = real(complex_fourier_coeffs(n));
                elseif strcmp(obj.boundary_condition, 'Closed')
                    sineCoeffs(n) = imag(complex_fourier_coeffs(n));
                elseif strcmp(obj.boundary_condition, 'Periodic')
                    cosineCoeffs(n) = real(complex_fourier_coeffs(n));
                    sineCoeffs(n) = imag(complex_fourier_coeffs(n));
                end
            end
        end

        %Function to reconstruct common phase in real space for a given
        %spectrum
        function out_com_phase = reconstruct_com_phase(obj, cosineCoeffs, sineCoeffs, z_grid)
            if nargin < 4
                z_grid = obj.z_grid;
            end

            ext_cutoff = min(length(cosineCoeffs), length(sineCoeffs));

            out_com_phase = zeros(1,length(z_grid));
            for i = 1:length(z_grid)
                com_phase = 0;
                for p = 1:ext_cutoff
                    com_phase = com_phase + cosineCoeffs(p)*cos(obj.wavevec_k*p*z_grid(i))+...
                        sineCoeffs(p)*sin(obj.wavevec_k*p*z_grid(i));
                end
                out_com_phase(i) = com_phase;
            end        
        end

        function out_com_phase = fdm_com_phase(obj, phase_init)
            dz = abs(obj.z_grid(2) - obj.z_grid(1));
            prefactor = (2*obj.m/(obj.hbar*obj.expansion_time))*(dz)^2;
            out_com_phase = zeros(1,length(obj.z_grid));
            out_com_phase(1) = phase_init(1);
            out_com_phase(2) = phase_init(2);
            for i = 3:length(obj.z_grid)
                disp(i)
                out_com_phase(i) = 2*out_com_phase(i)-out_com_phase(i-1) + prefactor*obj.dimensionless_density_ripple(i);
            end
        end
    end %End method
    methods (Static)
        %coherence factor fidelity
        function fcoh = fidelity_coh(phase_profile_1, phase_profile_2)
            fcoh = 0;
            phase_residue = phase_profile_1 - phase_profile_2;
            for i = 1:length(phase_profile_1)
                fcoh = fcoh + exp(1j*phase_residue(i));
            end
            fcoh = abs(fcoh/length(phase_residue));
        end
        %dot fidelity
        function fdot = fidelity_dot(phase_profile_1, phase_profile_2)
            fdot = dot(phase_profile_1,phase_profile_2)/...
                (norm(phase_profile_1).*norm(phase_profile_2));
        end

        %Squared error
        function aerr = mean_squared_error(phase_profile_1, phase_profile_2)
            aerr = mean((phase_profile_1 - phase_profile_2).^2);  
        end
     end
end %End class