classdef class_physical_parameters < handle
   
    properties
       % fundamental parameters
        c = 3e8; %speed of light
        m = (86.909*931.494e6)/(3e8*3e8); %mass of Rubidium atoms (eV s^2/m^2)
        hbar = 6.5821e-16; %reduced Planck constant (eV s)
        kb = 8.6173*1e-5; %Boltzmann constant in eV K^-1;

        %setting up default condensate parameters, geometry, and tof
        max_longitudinal_density = 75e6; %peak longitudinal density 75 atoms per microns
        omega = 2*pi*2e3; %x oscillator frequencies - transverse (s^-1)
        %omega = 3e3; 
        scattering_length = 5.2e-9; %scattering length for interaction broadening
        
        %by default, we have a condensate of length 100 microns, separated
        %a distance 3 microns initially, and spread over x_{max} - x_{min} microns after
        %tof
        default_condensate_length = 100e-6;%100 microns
        x_min = -50e-6; 
        x_max = 50e-6;
        default_separation_distance = 3e-6; %3 microns
        default_expansion_time = 15e-3; %15 ms
        default_coupling_J = 2*pi; %default tunnel coupling 2*pi*1Hz
        %default_coupling_J = 1;

        %default insitu density profile string 
        default_insitu_density = 'InverseParabola';
        
    end
    
end