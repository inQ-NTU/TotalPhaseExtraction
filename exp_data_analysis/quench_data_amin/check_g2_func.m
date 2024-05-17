clear all
close all
load('uncoupled_to_coupled_DW_box.mat')


%-2 ms data (beginning of the cooling), uncoupled double well (J = )
density_ripple = density_arr{:,2}; 
mean_density = mean_density_per_um{:,2};
phase_data = phase_arr{2,:};

%Compute g2 function
num_pix = size(density_ripple,2);
num_samples = size(density_ripple, 1);


g2 = zeros(num_pix, num_pix);
for i = 1:num_pix 
    for j = 1:num_pix
        for k = 1:num_samples
            g2(i,j) = g2(i,j)+density_ripple(k,i)*density_ripple(k,j)/(mean_density(i)*mean_density(j));
        end
    end
end
g2 = g2/num_samples;
