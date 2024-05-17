function psi_TOF = func_free_TOF_propagation_1D(psi,grid_si,TOF_si) 

%% constants
hbar     = 6.626e-34/(2*pi);      % Js
amuKg    = 1.6605387e-27;         % kg
Aatom = 87;                 % Rb atoms
massAtom = amuKg*Aatom; 

%% TOF propagation
no_gridpoints = length(grid_si);

dp = hbar*2.0*pi/(grid_si(end)-grid_si(1));
p = linspace(-floor(no_gridpoints/2)-1,floor(no_gridpoints/2)-1,no_gridpoints).*dp;
P=fftshift(p.^2/(2*massAtom));
Pop=exp(-1i*P*TOF_si./hbar);

psi_TOF=ifft(Pop.*fft(psi'));

end

