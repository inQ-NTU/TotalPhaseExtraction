classdef class_1d_correlation < handle
    properties
        %input
        all_phase_profiles

        %inferred
        referenced_phase_profiles
        average_phase_profile
        longitudinal_resolution
        nmb_of_sampled_profiles

        %output
        cov_matrix
    end
    methods 
        %implementing the constructor
        function obj = class_1d_correlation(phase_profile_matrix)
            obj.all_phase_profiles = phase_profile_matrix;
            obj.longitudinal_resolution = size(phase_profile_matrix,2);
            obj.nmb_of_sampled_profiles = size(phase_profile_matrix, 1);
            
            %calculating average phase profile - generally in our
            %simulation we set the average phase to be zero at every point,
            %but in this class, we do not assume that
            avg_phase_profile = zeros(1,obj.longitudinal_resolution);
            for i = 1:obj.nmb_of_sampled_profiles
                for j = 1:obj.longitudinal_resolution
                    avg_phase_profile(j) = avg_phase_profile(j)+phase_profile_matrix(i,j);           
                end
            end
            avg_phase_profile = (1/obj.nmb_of_sampled_profiles)*avg_phase_profile;
            obj.average_phase_profile = avg_phase_profile;
        end

        %Function to reference the phase profiles (set the midpoint to
        %zero)
        function ref_phases = reference_phase_profiles(obj, phase_profiles_data)
            if nargin < 2
                phase_profiles_data = obj.all_phase_profiles;
            end

            ref_phases = zeros(size(phase_profiles_data));
            for i = 1:obj.nmb_of_sampled_profiles
                ref_phases(i,:) = phase_profiles_data(i,:) - phase_profiles_data(i, floor(obj.longitudinal_resolution/2));
            end
            obj.referenced_phase_profiles = ref_phases; 
        end

        %Function to compute fourth order correlation G(z1, z2, z3, z4) function for fixed z3
        %and z4
        function corr = fourth_order_corr(obj, z3_idx, z4_idx, phase_profiles_data)
            if nargin<4
                phase_profiles_data = obj.all_phase_profiles;
            end

            %referencing the phase profiles
            for i = 1:obj.nmb_of_sampled_profiles
                phase_profiles_data(i,:) = phase_profiles_data(i,:) - phase_profiles_data(i,floor(obj.longitudinal_resolution/2)); 
            end

            %Computing the correlation
            corr = zeros(obj.longitudinal_resolution, obj.longitudinal_resolution);
            for i = 1:obj.longitudinal_resolution
                for j = 1:obj.longitudinal_resolution
                    prod = 0;
                    for k = 1:obj.nmb_of_sampled_profiles
                        prod = prod + phase_profiles_data(k,i)*phase_profiles_data(k,j)*phase_profiles_data(k, z3_idx)*phase_profiles_data(k, z4_idx);
                    end
                    prod = prod/obj.nmb_of_sampled_profiles;
                    corr(i,j) = prod;
                end
            end
        end
        
        %Function to compute correlation function of any order
        function corr = correlation_func(obj, order, phase_profiles_data)
            if nargin<3
                phase_profiles_data = obj.all_phase_profiles;
            end    
            %calculate all possible combination (and permutation)
            %of (z_1, z_2, ..., z_N) for N being the order
            combs = obj.nmultichoosek(1:obj.longitudinal_resolution, order);
            for i = 1:size(combs,1)
                unique_perms = unique(perms(combs(i,:)),'rows');
                combs = [combs;unique_perms(2:end,:)];
            end

            %Computing the average of the product of the phases for each
            %possible (z_1, z_2, ..., z_N)
            for i = 1:size(combs,1)
                avg = 0;
                for j = 1:obj.nmb_of_sampled_profiles
                    prod = 1;
                    for k = 1:order
                        prod = prod*(phase_profiles_data(j,combs(i,k))-phase_profiles_data(j, floor(obj.longitudinal_resolution/2)));
                    end
                    avg = avg+prod;
                end
                avg = avg/obj.nmb_of_sampled_profiles;
                idx = num2cell(combs(i,:));
                corr(idx{:}) = avg;
             end
        end
    
    %function to compute covariance matrix
    function cov_matrix = covariance_matrix(obj, phase_profiles_data)
        if nargin < 2
            phase_profiles_data = obj.all_phase_profiles;
        end
        %calculating covariance matrix
        cov_matrix = zeros(obj.longitudinal_resolution, obj.longitudinal_resolution);
        %averaging phi(z)phi(z') over all the samples
        for i = 1:obj.nmb_of_sampled_profiles
            for j = 1:obj.longitudinal_resolution
                for k =1:obj.longitudinal_resolution
                        cov_matrix(j,k) = cov_matrix(j,k)+(phase_profiles_data(i,j)-obj.average_phase_profile(j))*(phase_profiles_data(i,k)-obj.average_phase_profile(k));
                end
            end
        end
        cov_matrix = (1/obj.nmb_of_sampled_profiles)*cov_matrix;
        obj.cov_matrix = cov_matrix;
    end

    function w4 = wick_four_point_correlation(obj, phase_profiles_data)
        if nargin < 2
            phase_profiles_data = obj.all_phase_profiles;
        end
        %Compute two-point correlation function
        g2 = obj.correlation_func(2, phase_profiles_data);
    
        %Initialize four-point correlation function as a 4-tensor
        w4 = zeros(obj.longitudinal_resolution, obj.longitudinal_resolution, obj.longitudinal_resolution, obj.longitudinal_resolution);
        
        %Loop begin - compute the wick decomposition
        for i = 1:obj.longitudinal_resolution
            for j = 1:obj.longitudinal_resolution
                for k = 1:obj.longitudinal_resolution
                    for l = 1:obj.longitudinal_resolution
                        w4(i,j,k,l) = g2(i,j)*g2(k,l)+g2(i,k)*g2(j,l)+g2(i,l)*g2(j,k);
                    end
                end
            end
        end %Loop end
    end

    %function to compute g1 function
    function g1 = g1_corr(obj, phase_profiles_data)
        if nargin < 2
            phase_profiles_data = obj.all_phase_profiles;
        end
        g1 = zeros(obj.longitudinal_resolution, obj.longitudinal_resolution);
        for i = 1:obj.nmb_of_sampled_profiles
            for j = 1:obj.longitudinal_resolution
                for k = 1:obj.longitudinal_resolution
                    g1(j,k) = g1(j,k)+(cos(phase_profiles_data(i,j) - phase_profiles_data(i,k)));
                end
            end
        end
        g1  = (1/obj.nmb_of_sampled_profiles)*g1;
    end

    %function to compute correlation in Fourier space
    function fourier_cov = fourier_correlation(obj, phase_profiles_data)
        if nargin<2
            phase_profiles_data = obj.all_phase_profiles;
        end 
        %Compute fft for each profile
        fourier_data = transpose(abs(fft(transpose(phase_profiles_data)))*(1/size(phase_profiles_data,2)));
        
        %single sampling - remove frequency above Nyquist frequency
        l = floor(size(fourier_data,2)/2+1);
        fourier_data = fourier_data(:,1:l);

        %Compute correlation
        fourier_cov = zeros(l, l);
        for i = 1:obj.nmb_of_sampled_profiles
            for j = 1:l
                for k =1:l
                        fourier_cov(j,k) = fourier_cov(j,k)+fourier_data(i,j)*fourier_data(i,k);
                end
            end
        end
        fourier_cov = (1/obj.nmb_of_sampled_profiles)*fourier_cov;
        
    end
    end %end methods
  %%%%%%Static methods%%%%%%%%%%%%
  methods (Static)
        %function to compute combination (with repetition)
        %Code from https://stackoverflow.com/questions/28284671/generating-all-combinations-with-repetition-using-matlab
        function combs = nmultichoosek(values, k)
            if numel(values)==1
                n = values;
                combs = nchoosek(n+k-1,k);
            else
                n = numel(values);
                combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
                combs = reshape(values(combs),[],k);
            end
        end
        
        %function to compute distance between inferred and reference
        %covariance matrix
        %Two-norm (largest singular values)
        function dist = two_norm_distance(reference_cov, inferred_cov)
            dist = norm(reference_cov-inferred_cov)/norm(reference_cov);
        end

        %trace norm
        function dist = trace_norm_distance(reference_cov, inferred_cov)
            dist = sum(svd(reference_cov - inferred_cov))/sum(svd(reference_cov));
        end

        %Frobenius norm
        function dist = frobenius_norm_distance(reference_cov, inferred_cov)
            dist = norm(reference_cov-inferred_cov, 'fro')/norm(reference_cov, 'fro');
        end
  end %end static methods
end %end class
            