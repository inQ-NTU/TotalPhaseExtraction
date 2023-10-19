classdef class_relative_phase_extraction <  class_physical_parameters & handle

    properties

        %input
        separation_distance
        expansion_time
        input_tof_data

        %inferred
        longitudinal_resolution
        transversal_resolution
        normalized_tof_data

        %output
        phase_result_init
        phase_result_fit
        normalization_amplitudes
        contrasts
        gaussian_width
        all_fit_parameters
        reconstructed_interference_pattern

        %flag
        flag_fit_width

    end

    methods

        %% This implements the constructor
        function obj = class_relative_phase_extraction(input_tof_data, expansion_time, flag_fit_width, separation_distance)
            if nargin < 2
                obj.expansion_time = obj.default_expansion_time;
            else
                obj.expansion_time = expansion_time;
            end

            if nargin<3
                obj.flag_fit_width = 0;
            else
                obj.flag_fit_width = flag_fit_width;
            end

            if nargin<4
                obj.separation_distance = obj.default_separation_distance;
            else
                obj.separation_distance = separation_distance;
            end
            % 1. Store input data
            obj.input_tof_data = input_tof_data;
            obj.normalized_tof_data = obj.input_tof_data./max(obj.input_tof_data, [],'all');
            obj.longitudinal_resolution = size(input_tof_data, 1);
            obj.transversal_resolution = size(input_tof_data, 2);
        end

        %implement fitting formula independent of the interference class
        function rho_tof = transversal_expansion_formula(obj, x, amplitude, sigma, contrast, phase, shift)
            %wave number associated with fringe spacing
            kF = obj.separation_distance*(1+(obj.omega*obj.expansion_time)^2)/((sigma^2)*obj.omega*obj.expansion_time); 
            rho_tof = amplitude*exp(-x.^2/(sigma^2)).*(1+contrast*cos(kF*x+phase))+shift;       
        end

        %initial guess for fitting based on analytical formula in the
        %manuscript
        function guess = init_phase_guess(obj,rho_tof_data)

            if nargin < 2
                rho_tof_data = obj.normalized_tof_data;
            end

            calibration_slope = -obj.m*obj.separation_distance/(obj.hbar*obj.expansion_time);
            phase_guess = @(x) mod(calibration_slope*x,2*pi);
            slice_nmb = size(rho_tof_data, 1);
            pixel_nmb_x = size(rho_tof_data,2);
            guess = zeros(1,slice_nmb);

            for i = 1:slice_nmb
                %do the offsetting here
                input_interference_slice = rho_tof_data(i,:);
                [~, idx_max] = max(input_interference_slice);
                slice_maxima = obj.x_min+(idx_max/pixel_nmb_x)*(obj.x_max-obj.x_min);
                guessed_phase = phase_guess(slice_maxima);
                if guessed_phase > pi
                    guessed_phase = guessed_phase - 2*pi;
                end
                guess(i) = guessed_phase;
            end
            %Unwrapping process
            guess= obj.clean_phase_jump(guess);
            obj.phase_result_init = guess;
        end

        %fitting based on the transversal tof formula
        function fitted_phase = fitting(obj, guessed_phase_profile, rho_tof_data)

            if nargin < 3
                rho_tof_data = obj.normalized_tof_data;
            end
            %Define an objective/cost function from the interference class
            % write a loop for slice
            
            fitted_phase = zeros(1,obj.longitudinal_resolution);
            amp = zeros(1,obj.longitudinal_resolution);
            contr = zeros(1,obj.longitudinal_resolution);
            fit_params = cell(1, obj.longitudinal_resolution);
            reconstruced_tof = zeros(obj.longitudinal_resolution, obj.transversal_resolution);
            
            grid_x = linspace(obj.x_min, obj.x_max, obj.transversal_resolution);
            sigma_t = sqrt(obj.hbar/(obj.m*obj.omega))*sqrt(1+(obj.omega*obj.expansion_time)^2)*1e6; %in microns
            gauss_width = sigma_t*ones(1,obj.longitudinal_resolution);

            % Fmincon constraints: umplitude, contrast, phase
            if obj.flag_fit_width == 0
                search_lower_bound = [0,0,-inf,-2*pi];
                search_upper_bound = [inf,2,inf,2*pi];
            else
                search_lower_bound = [0,0,-2*pi,-inf,0];
                search_upper_bound = [inf,2, 2*pi,inf,100];
            end
            %Do fitting for each slice of image (fixed z)
            for i=1:obj.longitudinal_resolution
                interference_slice = rho_tof_data(i,:);
                if obj.flag_fit_width == 0
                    fitted_interference_slice = @(p) obj.transversal_expansion_formula(grid_x, p(1), sigma_t*1e-6, p(2), p(3),p(4));
                    tof_basic_cost_func = @(p) norm(interference_slice - fitted_interference_slice(p));
                    options = optimset('Display','none');
                    init_guess = [1,1,guessed_phase_profile(i),0];        
                    output = fmincon(tof_basic_cost_func, init_guess,[],[],[],[],search_lower_bound, search_upper_bound,[],options);
                    fitted_phase(i) = output(3);
                    amp(i) = output(1);
                    contr(i) = output(2);
                    reconstruced_tof(i,:) = fitted_interference_slice(output);
                    fit_params{i} = output;
                else
                    fitted_interference_slice = @(p) obj.transversal_expansion_formula(grid_x, p(1), p(5)*1e-6, p(2), p(3),p(4));
                    tof_basic_cost_func = @(p) norm(interference_slice - fitted_interference_slice(p));
                    options = optimset('Display','none');
                    init_guess = [1,1,guessed_phase_profile(i),0,sigma_t];        
                    output = fmincon(tof_basic_cost_func, init_guess,[],[],[],[],search_lower_bound, search_upper_bound,[],options);
                    fitted_phase(i) = output(3);
                    amp(i) = output(1);
                    contr(i) = output(2);
                    gauss_width(i) = output(5);
                    reconstruced_tof(i,:) = fitted_interference_slice(output);
                    fit_params{i} = output;
                end
            end
            fitted_phase = obj.clean_phase_jump(fitted_phase);
            obj.phase_result_fit = fitted_phase;
            obj.normalization_amplitudes = amp*max(obj.input_tof_data, [],'all');
            obj.contrasts = contr;
            obj.gaussian_width = gauss_width;
            obj.reconstructed_interference_pattern = reconstruced_tof*max(obj.input_tof_data,[],'all');
            obj.all_fit_parameters = fit_params;
        end

    end

    methods (Static)

        %dot product fidelity
        function fdot = fidelity_dot(phase_profile_1, phase_profile_2)

            fdot = dot(phase_profile_1,phase_profile_2)/...
                (norm(phase_profile_1).*norm(phase_profile_2));

        end

        %coherence factor fidelity
        function fcoh = fidelity_coh(phase_profile_1, phase_profile_2)
            fcoh = 0;
            phase_residue = phase_profile_1 - phase_profile_2;
            for i = 1:length(phase_profile_1)
                fcoh = fcoh + exp(1j*phase_residue(i));
            end
            fcoh = abs(fcoh/length(phase_residue));
        end

        %z dependent cosine distance
        function fcos = fidelity_cos_z_dependent(phase_profile_1, phase_profile_2)
            phase_residue = phase_profile_1 - phase_profile_2;
            fcos = cos(phase_residue);
        end


        %Create a function for phase unwrapping subroutine
        function output_phase_profile = clean_phase_jump(input_phase_profile)
                %Splitting data into two
                ref_index = ceil(length(input_phase_profile)/2);
                data_1 = input_phase_profile(1:ref_index);
                data_2 = input_phase_profile(ref_index:length(input_phase_profile));
                %forward unwrapping
                data_2 = unwrap(data_2);
                %backward unwrapping
                data_1 = flip(unwrap(flip(data_1)));
                %remove the reference point
                data_1(end) = [];
                %Combine the two
                output_phase_profile = [data_1 data_2];
        end
    end
end