classdef class_common_phase_extraction < class_physical_parameters & handle
    properties
        %input
        density_ripple
        mean_density
        expansion_time
        z_grid
        boundary_cut_ratio


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
        function obj = class_common_phase_extraction(density_ripple, mean_density, z_grid, condensate_length, expansion_time)
            if nargin < 3
                obj.mean_density = 2*obj.max_longitudinal_density; 
            else
                obj.mean_density = mean_density;
            end

            if nargin < 4
                obj.expansion_time = obj.default_expansion_time;
            else
                obj.expansion_time = expansion_time;
            end

            obj.z_grid = z_grid;
            obj.condensate_length = condensate_length;
            obj.mean_density_gradient = gradient(obj.mean_density);
            obj.wavevec_k = 2*pi/obj.condensate_length;

            obj.density_ripple = density_ripple;
            obj.lt = sqrt(obj.hbar*obj.expansion_time/(2*obj.m)); 
            obj.dimensionless_grid = obj.z_grid/obj.lt;
            obj.dimensionless_mean_density_gradient = gradient(obj.mean_density, obj.dimensionless_grid)./obj.mean_density;
            
            dimensionless_density_ripple = 1-(obj.density_ripple./obj.mean_density);
            dz = obj.z_grid(2) - obj.z_grid(1);
            obj.dimensionless_density_ripple = dimensionless_density_ripple - (dz/obj.condensate_length)*sum(dimensionless_density_ripple);
            
        end
        

        %2. Function to extract common phase spectrum when the mean density is
        %uniform
        function [cosineCoeffs, sineCoeffs] = extract_com_spectrum_uniform(obj, ext_cutoff)

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
                cosineCoeffs(n) = real(complex_fourier_coeffs(n));
                sineCoeffs(n) = imag(complex_fourier_coeffs(n));
            end
        end

        %Function to interpolate ODE coefficients corresponding to
        %continuity equation
        function [A,B] = interp_ode_coefficients(obj, nu)
            A = interp1(obj.dimensionless_grid,obj.dimensionless_mean_density_gradient, nu);
            B = interp1(obj.dimensionless_grid, obj.dimensionless_density_ripple,  nu);
        end
        %Function to setup continuity equation - not used in the current
        %manuscript
        function dxdnu = continuity_eqn(obj, nu, x)
            [A,B] = obj.interp_ode_coefficients(nu);
            dxdnu =-A*x+B;
        end

        %Solving for velocity profile
        function [vel_profile, z_grid] = extract_common_velocity(obj, nugrid, boundary_condition)
            if nargin < 2
                nugrid= obj.dimensionless_grid;
            end

            if nargin < 3
                boundary_condition = 0;
            end 
            [nu,x] = ode45(@(nu,x) obj.continuity_eqn(nu, x), nugrid, boundary_condition);
            z_grid = nu*obj.lt;
            vel_profile = x/obj.lt;
        end

        %Solving for common phase spectrum from the obtained velocity profile 
        % assuming the absence of zero mode
        function [cosineCoeffs, sineCoeffs] = extract_com_spectrum_continuity(obj, ext_cutoff)
            [vel_profile, intermediate_z_grid] = obj.extract_common_velocity();
            
            
            %Eliminate zero mode in the velocity profile
            dz = intermediate_z_grid(2)-intermediate_z_grid(1);
            vel_profile = vel_profile - (dz/obj.condensate_length)*sum(vel_profile);

            cosineCoeffs = zeros(1, ext_cutoff);
            sineCoeffs = zeros(1,ext_cutoff);
            %Perform Fourier analysis of the common velocity
            for p = 1:ext_cutoff
                fp = 0;
                for i = 1:length(intermediate_z_grid)
                    fp = fp + vel_profile(i)*exp(-1j*p*obj.wavevec_k*intermediate_z_grid(i));
                end
                fp = fp*(dz/(pi*p));
                cosineCoeffs(p) = imag(fp);
                sineCoeffs(p) = real(fp);
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