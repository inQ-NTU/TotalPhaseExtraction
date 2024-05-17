function g0 = calc_g0(omegap_si)

hqsi = 1.054571726e-34; % reduced Planck's constant in SI-units
as_si = 5.2e-9; % scattering length in m

% calculate the interaction constant
g0 = 2*hqsi*as_si*omegap_si;

end