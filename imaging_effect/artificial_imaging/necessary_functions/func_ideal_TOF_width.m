function width = func_ideal_TOF_width(TOF_si,omega_trans)

% constants
hbar     = 6.626e-34/(2*pi);      % Js
amuKg    = 1.6605387e-27;         % kg
Aatom = 87;                 % Rb atoms
massAtom = amuKg*Aatom;

% free expansion transversal
width0=sqrt(hbar/(massAtom*omega_trans)); %in trap width = harmonic oscillator ground state width
width=width0*sqrt(1+(TOF_si*omega_trans)^2); %after free expansion

end