classdef class_interference_pattern< class_physical_parameters

    properties

        %input phase profile
        phase_profile_RS
        relative_phase_profile
        common_phase_profile
        get_relative_phase_z
        get_common_phase_z

        %tunable parameters
        expansion_time
        separation_distance_d
        condensate_length_Lz
        buffer_length

        %flags
        flag_interaction_broadening  %0: no interaction broadening, 
        % 1: include interaction broadening
        
        %insitu density
        insitu_density_profile_str
        insitu_density

        %grids - there are three grids:
        %1. input grid - determined by input phase profile
        %2. output grid - taking into account both input and buffer (for
        %longitudinal expansion)
        %3. integration grid - a fine grid to do integration
        nmb_input_points_z
        nmb_integration_points_z
        nmb_buffer_points_z
        nmb_output_points_z
        nmb_output_points_x

        input_grid_z
        input_grid_x
        integration_grid_z
        output_grid_z
        output_grid_x

        transversal_density_avg_squared
        phase_shift_d_x_t

    end


    methods

        %% Constructor
        function obj = class_interference_pattern(phase_profile_RS, expansion_time, ...
                insitu_density_profile_str, flag_interaction_broadening, ...
                separation_distance_d, condensate_length_Lz, buffer_length)

            % 1. Run the parent class constructor
            obj = obj@class_physical_parameters();
            % 2. Set up some of the physics parameters, possibly to default
            if nargin < 2
                obj.expansion_time = obj.default_expansion_time;
            else
                obj.expansion_time = expansion_time;
            end

            if nargin < 3
                obj.insitu_density_profile_str = obj.default_insitu_density;
            else
               obj.insitu_density_profile_str = insitu_density_profile_str;     
            end

            if nargin < 4
                obj.flag_interaction_broadening = 0;
            else
                obj.flag_interaction_broadening = flag_interaction_broadening;
            end

            if nargin < 5
                obj.separation_distance_d = obj.default_separation_distance;
            else
                obj.separation_distance_d = separation_distance_d;
            end

            if nargin <6
                obj.condensate_length_Lz = obj.default_condensate_length;
            else
                obj.condensate_length_Lz = condensate_length_Lz;
            end

            if nargin<6
                obj.buffer_length = 0;
            else
                obj.buffer_length = buffer_length;
            end

            % 4. Set up profile, if dimension is 2 -> relative + common
            %if dimension is 1 -> relative
            if size(phase_profile_RS,1) == 2
                obj.phase_profile_RS = phase_profile_RS;
                obj.relative_phase_profile = phase_profile_RS(1,:);
                obj.common_phase_profile = phase_profile_RS(2,:);
            elseif size(phase_profile_RS,1) == 1
                obj.phase_profile_RS = phase_profile_RS;
                obj.relative_phase_profile = phase_profile_RS;
                obj.common_phase_profile = zeros(size(phase_profile_RS));
            end

            % 5. Grids
            % 5.1 Input grid (matching the phase profile)
            obj.nmb_input_points_z = length(obj.relative_phase_profile);
            obj.input_grid_z =  linspace(-obj.condensate_length_Lz/2, obj.condensate_length_Lz/2, obj.nmb_input_points_z );

            % 5.2 Output grid
            obj.nmb_output_points_z = obj.nmb_input_points_z;
            obj.nmb_buffer_points_z = floor((obj.buffer_length/obj.condensate_length_Lz)*obj.nmb_output_points_z);
            obj.output_grid_z = linspace(-obj.buffer_length-obj.condensate_length_Lz/2,...
                (obj.condensate_length_Lz/2)+obj.buffer_length, obj.nmb_output_points_z+2*obj.nmb_buffer_points_z);
            
            %5.3 Integration grid
            obj.nmb_integration_points_z = floor((2.6/obj.expansion_time)+238.23); %phenomenological formula
            obj.integration_grid_z = linspace(-obj.condensate_length_Lz/2, obj.condensate_length_Lz/2, obj.nmb_integration_points_z);
            
            % 5.4 Setup of transversal grid
            
            % By default we have 160um transversally and 100um
            % longitudinally so we use proportionally to roughly get a
            % square grid
            obj.nmb_output_points_x = ceil( ...
                (obj.x_max-obj.x_min) / obj.condensate_length_Lz ...
                * obj.nmb_input_points_z );
            obj.output_grid_x = linspace(obj.x_min, obj.x_max,...
                obj.nmb_output_points_x);

            % 8. Prepare for integrals
            % 8.1 Interpolate phase profiles using condensate grid
            obj.get_relative_phase_z = obj.phase_profile_interpolation_init(obj.input_grid_z, obj.relative_phase_profile);
            obj.get_common_phase_z = obj.phase_profile_interpolation_init(obj.input_grid_z, obj.common_phase_profile);

            %9. Assign insitu density
            density = zeros(1, obj.nmb_input_points_z);
            for i = 1:obj.nmb_input_points_z
                density(i) = obj.longitudinal_density(obj.input_grid_z(i));
            end
            obj.insitu_density = density;
        
        end

        %Expansion dynamics
        %Initial width of the Gaussian - can be with or without interaction
        %broadening
        function density_sigma_init = compute_density_sigma_init(obj, z)
            if obj.flag_interaction_broadening == 0
                density_sigma_init = sqrt( obj.hbar/(obj.m*obj.omega) );
            else
                density_sigma_init = sqrt( obj.hbar/(obj.m*obj.omega) )*sqrt(1+obj.scattering_length*obj.longitudinal_density(z));
            end
        end

        function density_sigma_t = compute_density_sigma_t(obj, z, time)
            if nargin <3
                time = obj.expansion_time;
            end
            density_sigma_init = obj.compute_density_sigma_init(z);
            density_sigma_t = density_sigma_init*sqrt(1+obj.omega^2*time.^2);
       end
        
        % single particle Green's function
        function single_particle_green_function = green_function(obj,z,zp,t)
            single_particle_green_function = sqrt( obj.m/(2*pi*1j*obj.hbar*t) )...
                *exp( ( (obj.m)/(2*1j*obj.hbar*t) )*(z-zp).^2 );
        end

        %% TOF interference functions

        % phase prefactors used in 3D ToF expansion -> reduce the number of
        % necessary function
        function transversal_density = integrand_density_prefactor(obj,separation_distance_d,longitudinal_position_z,expansion_time)
            if nargin < 4
                expansion_time = obj.expansion_time;
            end
            density_sigma_t = obj.compute_density_sigma_t(longitudinal_position_z, expansion_time);
            transversal_density = class_interference_pattern.normalized_Gaussian(obj.output_grid_x,...
                separation_distance_d/(2), density_sigma_t);
        end

        function fringe_spacing = compute_fringe_spacing(obj, separation_distance_d, z)
            fringe_spacing = 2*pi*(obj.compute_density_sigma_init(z)^2)*obj.omega*obj.expansion_time/separation_distance_d;
        end

        function tof_dependent_phase_factor = x_dependent_integrand_prefactor(obj, separation_distance_d, z)
            tof_dependent_phase_factor = exp(((1j*pi)/(obj.compute_fringe_spacing(separation_distance_d, z))).*obj.output_grid_x );
        end

        % phase prefactors used in transversal expansion
        function transversal_density_avg_squared = compute_density_prefactor(obj, z)
            density_prefactor_right = obj.integrand_density_prefactor(obj.separation_distance_d,z,obj.expansion_time);
            density_prefactor_left = obj.integrand_density_prefactor(-obj.separation_distance_d,z,obj.expansion_time);
            transversal_density_avg_squared = (1/2)*(density_prefactor_left.^2 + density_prefactor_right.^2);
        end

        function phase_shift_d_x_t = compute_phase_shift_d_x_t(obj,separation_distance_d, z)
            phase_shift_d_x_t = ( (2*pi)/obj.compute_fringe_spacing(separation_distance_d, z))*(obj.output_grid_x);
        end

        %Longitudinal Density profiles
        %Thomas-Fermi density - Inverse parabola
        function rho = inverse_parabola_density(obj, z)
            rho = obj.max_longitudinal_density(1-(2*z/obj.condensate_length_Lz)^2)*obj.step_func((1/2)-abs(z));
        end

        %Thomas-Fermi density - box potential
        function rho = box_density(obj, z)
            gas_length = obj.condensate_length_Lz; %Gas length in microns
            max_density = obj.max_longitudinal_density;
            falloff_param = 0.9;
            rho = (max_density/2)*(tanh(z*1e6+falloff_param*gas_length*1e6/2)-tanh(z*1e6-falloff_param*gas_length*1e6/2));
        end

      
        %In case if we want flat density profile
        function rho = flat_density(obj, z)
            rho = obj.max_longitudinal_density;
        end

        %Final longitudinal density
        function rho = longitudinal_density(obj, z)
            if strcmp(obj.insitu_density_profile_str, 'InverseParabola')
                rho = obj.inverse_parabola_density(z);
            elseif strcmp(obj.insitu_density_profile_str, 'BoxPotential')
                rho = obj.box_density(z);
            elseif strcmp(obj.insitu_density_profile_str, 'Flat')
                rho = obj.flat_density(z);
            else
                rho = obj.inverse_parabola_density(z);
                disp('The insitu density is not known. Using the default inverse parabola instead.')
            end
        end

        %TOF FORMULA: transversal expansion only%
        function rho_tof = tof_transversal_expansion(obj, relative_phase)
            if nargin < 2
                relative_phase = obj.get_relative_phase_profile;
            end
            rho_tof = zeros(length(obj.output_grid_z), length(obj.output_grid_x));
            for j = obj.nmb_buffer_points_z+1:length(obj.input_grid_z)+obj.nmb_buffer_points_z
                rho_tof(j,:) = 2*obj.longitudinal_density(obj.output_grid_z(j))*obj.compute_density_prefactor(obj.output_grid_z(j)).*...
                    ( 1 + cos(relative_phase(j-obj.nmb_buffer_points_z) ...
                + obj.compute_phase_shift_d_x_t(obj.separation_distance_d, obj.output_grid_z(j))) );
            end
        end

        %TOF FORMULA - full with longitudinal expansion (using Riemann sum
        %for integration
        function rho_tof = tof_full_expansion(obj, relative_phase_profile, common_phase_profile)
            if nargin < 3
                get_relative_phase_z = obj.get_relative_phase_z;
                get_common_phase_z = obj.get_common_phase_z;
            end
            if nargin == 3
                get_relative_phase_z = obj.phase_profile_interpolation_init(obj.input_grid_z, ...
                    relative_phase_profile);
                get_common_phase_z = obj.get_common_phase_z;
            end

            integrand_for_any_z = @(z_prime)  (...
                ( obj.green_function( obj.output_grid_z, z_prime, obj.expansion_time )') .* ...
                ( obj.integrand_density_prefactor( obj.separation_distance_d,z_prime).* ...
                obj.x_dependent_integrand_prefactor(obj.separation_distance_d, z_prime)...
                * sqrt(obj.longitudinal_density(z_prime))...
                *exp( (1j/2) * get_common_phase_z( z_prime ) ) ...
                * exp( (1j/2) *get_relative_phase_z( z_prime ) ) ...
                + ...
                obj.integrand_density_prefactor( -obj.separation_distance_d, z_prime) .*...
                obj.x_dependent_integrand_prefactor( - obj.separation_distance_d, z_prime)...
                *sqrt(obj.longitudinal_density(z_prime))...
                * exp( (1j/2)*get_common_phase_z(z_prime) ) ...
                *  exp( -(1j/2)* get_relative_phase_z( z_prime) )...
                ));

            %Rectangle Method:
            grid_spacing_integral = obj.integration_grid_z(2) - obj.integration_grid_z(1);
            integral = 0 ;
            for z_prime = obj.integration_grid_z
                integral = integral + integrand_for_any_z(z_prime)*grid_spacing_integral;
            end
            rho_tof = abs(integral).^2;
        end


       %Processing of the TOF interference image

      %Normalizing the interference pattern - \int rho dx dz = N where N is
      %the total number of atoms
      function normalized_rho_tof = normalize(obj, rho_tof, N)
          sum_density = sum(rho_tof, "all");
          pixel_width_z = obj.condensate_length_Lz/size(rho_tof, 1);
          pixel_width_x = (obj.x_max - obj.x_min)/size(rho_tof, 2);
          normalizing_coeff = sum_density*pixel_width_x*pixel_width_z;
          normalized_rho_tof = (N/normalizing_coeff)*rho_tof;
      end


        function convolved_tof = convolution2d(obj, rho_tof, convolution_scale)
            if nargin <3
                convolution_scale = obj.default_2d_conv_scale;
            end
            [m,n] = size(rho_tof);
            %calculating sigma relative to the system size
            sigma_m = convolution_scale * m;
            sigma_n = convolution_scale * n;

            S_m = class_interference_pattern.convolution_matrix(sigma_m, m);
            S_n = class_interference_pattern.convolution_matrix(sigma_n, n);

            convolved_tof = S_m*rho_tof*S_n;

        end

        % get interpolated phase profile
        function relative_phase_profile = get_relative_phase_profile(obj)
            relative_phase_profile = obj.get_relative_phase_z( obj.input_grid_z );
        end
    end

    methods (Static)

        function phase_profile_interpolant = phase_profile_interpolation_init(z_grid, phase)
            [Z] = ndgrid(z_grid);
            phase_profile_interpolant = griddedInterpolant( Z, phase );
        end

      function S = convolution_matrix( sigma, nmb_grid_points)
            % Calculates the convolution matrix
            S = zeros( nmb_grid_points );
            for x = 1:nmb_grid_points
                S( x, : ) = exp( - (  (1:nmb_grid_points) - x ).^2 / (2 * sigma^2) );
            end
            %normalising the convolution matrix
            for x = 1:nmb_grid_points
                norm = sum( S(x,:));
                S( x, : ) = S( x,: ) / norm;
            end
      end
      
      function f = normalized_Gaussian(x, mean, standard_deviation)
            %the normalized Gaussian functions are given as
            f = (1/(pi*standard_deviation^2)^(1/4))*exp(- (x - mean).^2/(2*standard_deviation^2)) ;
      end

      %Step function for defining longitudinal density
      function theta = step_func(x)
        if x>0
            theta = 1;
        elseif x<0
            theta = 0;
        else
            theta = 0.5;
        end
      end

      %function to add Gaussian noise + baseline noise
      function noisy_rho_tof = add_gaussian_noise(rho_tof_data, sigma_scale)
            max_tof_value = max(rho_tof_data, [], 'all');
            sigma = sigma_scale*max_tof_value; %the std is the scaling*absolute peak of the data
            noise_data = randn(size(rho_tof_data))*sigma; %generate gaussian random noise
            noisy_rho_tof = rho_tof_data+noise_data; %additive noise

            %set a baseline value to prevent negative intensity
            baseline = min(noisy_rho_tof,[],'all');
            if baseline<0
                noisy_rho_tof = noisy_rho_tof-baseline;
            end
      end

    end
end