function [phase_distribution, dummy]=DensityDependantOU_SingleQuasiBECFun(T, ~, dz, n1d)

% This function uses an updating formula for the Ornstein-Uhlenbeck
% to process to simulate an in situ phase distribution of a single quasi
% condensate

%(see H.P. Stimming et al, PRL 105 015301 and reference 25 therein PRE 54, 2084)
% In every spae step the parameters of the stochastic process can be
% changed. This allows e.g to account for varying density!

% Vectorized version.

% TS 2016-08 corrected a small bug, until now always startphase 0 was used for the random process 


%% constants
hbar     = 6.626e-34/(2*pi);      % Js
amuKg    = 1.6605387e-27;         % kg
kB       = 1.381e-23;             % J/K
Aatom = 87;                 % Rb atoms
massAtom = amuKg*Aatom;
dummy = 0;                  % for downward compatibility, since the old OU code gives out a temporary profile no one uses

%% OU
kappaT = massAtom*kB*T/hbar^2./n1d;    % However, n1d = 0 at a number of grid points.
singularities = (kappaT == Inf);       % Find these gridpoints.
kappaT(singularities) = 0;             % Fix them. Only Chuck Norris can divide by zero.
phase_distribution = (cumsum(sqrt(kappaT*dz).*randn(1,length(n1d))))'; % Conjugated for downward compatibility.

%add a global random phase as the start value of the random process was
%always 0
phase_distribution = 2*pi*rand(1) + phase_distribution;

phase_distribution(singularities) = 0; % Without this, cumsum creates finite values for the phase for empty grid points at the "right" side of the cloud. This sets them to zero, as in the old script.

end
