classdef class_common_phase_spectrum < class_physical_parameters & handle
    properties
        %input
        density_ripple
        expansion_time
        z_grid
        condensate_length

        %intermediate variables
        wavevec_k
        
        %output
        com_phase_cosine_spectrum
        com_phase_sine_spectrum
        ripple_fourier_coeffs
        com_phase_profile
    end

    methods
        function obj = class_common_phase_spectrum(density_ripple_data, z_grid, expansion_time, condensate_length)
            if nargin < 3
                obj.expansion_time = obj.default_expansion_time;
            else
                obj.expansion_time = expansion_time;
            end

            if nargin < 4
                obj.condensate_length = obj.default_condensate_length;
            else
                obj.condensate_length = condensate_length;
            end

            %Store input data
            obj.density_ripple = density_ripple_data;
            obj.z_grid = z_grid;
            obj.wavevec_k = 2*pi/obj.condensate_length;
        end
        
        %Flag deconvolution -> mitigate the effect of rectangular windowing
        %Default value is 1 (activated), to deactivate simply set it to 0
        function [cosineCoeffs, sineCoeffs] = extract_com_spectrum(obj, n_max_fourier)

            %initialize the outputs
            cosineCoeffs = zeros(1, n_max_fourier);
            sineCoeffs = zeros(1,n_max_fourier);
            complex_fourier_coeffs = zeros(1,n_max_fourier);

            %Compute and initialize parameters
            lt = sqrt(obj.hbar*obj.expansion_time/obj.m);
            eps_t = lt/obj.condensate_length; 
            dz = obj.z_grid(2) - obj.z_grid(1);
            ripple = obj.density_ripple;
            
            %Making sure zero constant shift in Fourier decomposition
            ripple = ripple - sum(ripple)*dz/obj.condensate_length; 
            
            %Start computing Fourier coefficients
            for n = 1:n_max_fourier
                fn = 0;
                for i = 1:length(obj.z_grid)
                    fn = fn + ripple(i)*exp(-1j*n*obj.wavevec_k*obj.z_grid(i));
                end
                complex_fourier_coeffs(n) = (dz/obj.condensate_length)*fn;
            end

            %Assigning outputs
            obj.ripple_fourier_coeffs = complex_fourier_coeffs;
            for n = 1:n_max_fourier
                complex_fourier_coeffs(n) = -(1/((pi*eps_t*n)^2))*complex_fourier_coeffs(n);
                cosineCoeffs(n) = real(complex_fourier_coeffs(n));
                sineCoeffs(n) = -imag(complex_fourier_coeffs(n));
            end
            obj.com_phase_cosine_spectrum = cosineCoeffs;
            obj.com_phase_sine_spectrum = sineCoeffs;
        end

        %Reconstructing common phase in real space
        function com_phase_profile = extract_com_profile(obj, z_grid)
            com_phase_profile = zeros(1,length(z_grid));
            n_max_fourier = length(obj.com_phase_sine_spectrum);
            for i = 1:length(z_grid)
                for n = 1:n_max_fourier
                    cosine_term = obj.com_phase_cosine_spectrum(n)*cos(n*obj.wavevec_k*z_grid(i));
                    sine_term = obj.com_phase_sine_spectrum(n)*sin(n*obj.wavevec_k*z_grid(i));
                    com_phase_profile(i) = com_phase_profile(i)+cosine_term+sine_term;
                end
            end
            obj.com_phase_profile = com_phase_profile;
        end
    end
    
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
     end

end