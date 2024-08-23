%This is a code to simulate interference pattern data acquisition after
%time of flight expansion of parallel 1D Bose gases. The TOF dynamics is
%assumed to be fully ballistic (no final state interaction)
classdef class_interference_pattern < class_physical_parameters

    properties

        %input phase profile
        phase_profile_1
        phase_profile_2

        get_phase_profile_1
        get_phase_profile_2

        %tunable parameters
        expansion_time
        separation_distance_d
        condensate_length_Lz
        buffer_length_x
        buffer_length_z
        integration_buffer_length

        %flags
        flag_interaction_broadening  
        %0: no interaction broadening, 
        %1: include interaction broadening
        
        %insitu density
        insitu_density_1
        insitu_density_2
        get_insitu_density_1
        get_insitu_density_2

        %grids - there are three grids:
        %1. input grid - determined by input phase profile
        %2. output grid - output grid for TOF picture 
        %3. integration grid - a fine grid to do integration
        nmb_input_points_z
        nmb_integration_points_z
        nmb_buffer_points_z
        nmb_buffer_points_x
        nmb_output_points_z
        nmb_output_points_x

        input_grid_z
        integration_grid_z
        output_grid_z
        output_grid_x

    end


    methods

        %% Constructor
        function obj = class_interference_pattern(phase_profile_all, density_profile_all, expansion_time, ...
                buffer_length, condensate_length_Lz, flag_interaction_broadening, ...
                separation_distance_d)

            % 1. Run the parent class constructor
            obj = obj@class_physical_parameters();

            % 2. Set up some of the physics parameters, possibly to default
            if nargin < 3
                obj.expansion_time = obj.default_expansion_time;
            else
                obj.expansion_time = expansion_time;
            end

            if nargin<4
                obj.buffer_length_z = 0;
                obj.buffer_length_x = 0;
            else
                obj.buffer_length_z = buffer_length(1);
                obj.buffer_length_x = buffer_length(2);
            end

            if nargin <5
                obj.condensate_length_Lz = obj.default_condensate_length;
            else
                obj.condensate_length_Lz = condensate_length_Lz;
            end

            if nargin < 6
                obj.flag_interaction_broadening = 0;
            else
                obj.flag_interaction_broadening = flag_interaction_broadening;
            end

            if nargin < 7
                obj.separation_distance_d = obj.default_separation_distance;
            else
                obj.separation_distance_d = separation_distance_d;
            end

            %3. Set up input phase and density profiles
            if size(phase_profile_all,1) == 2
               obj.phase_profile_1 = phase_profile_all(1,:);
               obj.phase_profile_2 = phase_profile_all(2,:);
            end
            %If dimension is 2 -> gas 1 and 2
            %If dimension is 1 -> gas 1 = gas 2
            if size(density_profile_all,1) == 2
                obj.insitu_density_1 = density_profile_all(1,:);
                obj.insitu_density_2 = density_profile_all(2,:);
            elseif size(density_profile_all,1) == 1
                obj.insitu_density_1 = density_profile_all;
                obj.insitu_density_2 = density_profile_all;
            end

            % 4. Grids
            % 4.1 Input grid (matching the phase profile)
            obj.nmb_input_points_z = length(obj.phase_profile_1);
            obj.input_grid_z =  linspace(-obj.condensate_length_Lz/2, obj.condensate_length_Lz/2, obj.nmb_input_points_z );

            % 4.2 Output grid
            obj.nmb_output_points_z = obj.nmb_input_points_z;
            obj.nmb_buffer_points_z = floor((obj.buffer_length_z/obj.condensate_length_Lz)*obj.nmb_output_points_z);
            obj.output_grid_z = linspace(-obj.buffer_length_z-obj.condensate_length_Lz/2,...
                (obj.condensate_length_Lz/2)+obj.buffer_length_z, obj.nmb_output_points_z+2*obj.nmb_buffer_points_z);
            %Without buffer output_grid_x is assumed to be the same as
            %output_grid_z
            obj.nmb_buffer_points_x = floor((obj.buffer_length_x/obj.condensate_length_Lz)*obj.nmb_output_points_z);
            obj.output_grid_x = linspace(-obj.buffer_length_x-obj.condensate_length_Lz/2,...
                (obj.condensate_length_Lz/2)+obj.buffer_length_x, obj.nmb_output_points_z+2*obj.nmb_buffer_points_x);
            
            %4.3 Integration grid
            obj.integration_buffer_length = 0.1*obj.condensate_length_Lz+obj.buffer_length_z; %Default: integration buffer length 10% of the initial gas length
            obj.nmb_integration_points_z = 8*obj.nmb_input_points_z; %Default: Refine grid for integration, 8 times more refine than the input grid by default
            obj.integration_grid_z = linspace(-obj.condensate_length_Lz/2-obj.integration_buffer_length, ...
                obj.condensate_length_Lz/2+obj.integration_buffer_length, obj.nmb_integration_points_z);
            
            %4.4 Interpolating the profile
            obj.get_phase_profile_1 = obj.profile_interpolation_init(obj.input_grid_z, obj.phase_profile_1);
            obj.get_phase_profile_2 = obj.profile_interpolation_init(obj.input_grid_z, obj.phase_profile_2);
            obj.get_insitu_density_1 = obj.profile_interpolation_init(obj.input_grid_z, obj.insitu_density_1);
            obj.get_insitu_density_2 = obj.profile_interpolation_init(obj.input_grid_z, obj.insitu_density_2);
        end
        
        %5. Expansion dynamics
        %5.1 Initial width of the Gaussian - can be with or without interaction broadening
        function density_sigma_init = compute_density_sigma_init(obj, z)
            if obj.flag_interaction_broadening == 0
                density_sigma_init = sqrt(obj.hbar/(obj.m*obj.omega) );
            elseif obj.flag_interaction_broadening == 1
                density_sigma_init = sqrt(obj.hbar/(obj.m*obj.omega) )*sqrt(1+obj.scattering_length*obj.longitudinal_density(z));
            end
        end

        %5.2 Cloud width after expansion time t
        function density_sigma_t = compute_density_sigma_t(obj, z, time)
            if nargin <3
                time = obj.expansion_time;
            end
            density_sigma_init = obj.compute_density_sigma_init(z);
            density_sigma_t = density_sigma_init*sqrt(1+obj.omega^2*time.^2);
        end
        
        %5.3 Single particle Green's function
        function single_particle_green_function = green_function(obj,z,zp,t)
            single_particle_green_function = sqrt(obj.m/(2j*pi*obj.hbar*t))...
                *exp((1j*obj.m/(2*obj.hbar*t))*(z-zp).^2); 
        end

        %%6. TOF interference functions
        %6.1 Transversal field prefactor
        function transversal_field = transversal_field_prefactor(obj,separation_distance_d,longitudinal_position_z,...
                expansion_time)
            if nargin < 4
                expansion_time = obj.expansion_time;
            end
            density_sigma_t = obj.compute_density_sigma_t(longitudinal_position_z, expansion_time);
            transversal_field = obj.sqrt_normalized_Gaussian(obj.output_grid_x,...
                separation_distance_d/2, density_sigma_t);
        end

        %6.2 Fringe spacing lambda = 2 pi/k where k is the associated wavenumber
        function fringe_spacing = compute_fringe_spacing(obj, separation_distance_d, z)
            fringe_spacing = 2*pi*(obj.compute_density_sigma_init(z)^2)*obj.omega*obj.expansion_time/separation_distance_d;
        end

        %6.3 phase shift due to transverse expansion = kx/2 =pi x/lambda
        function alpha = tof_dependent_phase_prefactor(obj, separation_distance_d, z)
            alpha = exp((1j*pi/obj.compute_fringe_spacing(separation_distance_d, z)).*obj.output_grid_x );
        end

        %6.4: TOF EXPANSION FORMULA
        %Function for simulating TOF without longitudinal expansion
        %(transverse expansion only)
        function rho_tof = tof_transverse_expansion(obj)
            field_1 = @(z) obj.transversal_field_prefactor(obj.separation_distance_d, z).*obj.tof_dependent_phase_prefactor(obj.separation_distance_d, z)....
                *sqrt(obj.get_insitu_density_1(z))*exp(1j*obj.get_phase_profile_1(z));
            field_2 = @(z) obj.transversal_field_prefactor(-obj.separation_distance_d, z).*obj.tof_dependent_phase_prefactor(-obj.separation_distance_d, z)....
                *sqrt(obj.get_insitu_density_2(z))*exp(1j*obj.get_phase_profile_2(z));
            for i = 1:length(obj.output_grid_z)
                rho_tof(i,:) = abs(field_1(obj.output_grid_z(i)) + field_2(obj.output_grid_z(i))).^2;
            end
        end

        %Function for simulating TOF with longitudinal expansion
        function rho_tof = tof_full_expansion(obj)
            integrand_for_any_z = @(z_prime)transpose(obj.green_function(obj.output_grid_z, z_prime, obj.expansion_time)).* ...
                ( ...
                obj.transversal_field_prefactor(obj.separation_distance_d, z_prime).* ...
                obj.tof_dependent_phase_prefactor(obj.separation_distance_d, z_prime)...
                *sqrt(obj.get_insitu_density_1(z_prime))...
                *exp(1j*obj.get_phase_profile_1(z_prime))...
                + ...
                obj.transversal_field_prefactor(-obj.separation_distance_d, z_prime) .*...
                obj.tof_dependent_phase_prefactor(-obj.separation_distance_d, z_prime)...
                *sqrt(obj.get_insitu_density_2(z_prime))...
                *exp(1j*obj.get_phase_profile_2(z_prime))...
                );

            %Rectangle Method:
            grid_spacing_integral = abs(obj.integration_grid_z(2) - obj.integration_grid_z(1));
            integral = 0 ;
            for z_prime = obj.integration_grid_z
                integral = integral+integrand_for_any_z(z_prime)*grid_spacing_integral;
            end
            rho_tof = abs(integral).^2;
        end
  end
    methods (Static)

        function profile_interpolant = profile_interpolation_init(z_grid, phase)
            [Z] = ndgrid(z_grid);
            profile_interpolant = griddedInterpolant( Z, phase );
        end
      
        function f = sqrt_normalized_Gaussian(x, mean, standard_deviation)
            %the root of normalized Gaussian function is given as
            f = (1/(pi*standard_deviation^2)^(1/4))*exp(- (x - mean).^2/(2*standard_deviation^2)) ;
      end
  end
end