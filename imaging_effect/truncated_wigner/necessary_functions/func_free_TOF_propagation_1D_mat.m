function psi_TOF = func_free_TOF_propagation_1D_mat(psi,grid_si,TOF_si,space_dim) 

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

size_psi = size(psi);
no_dim_psi = length(size_psi);
reshape_vec = ones(1,no_dim_psi);
reshape_vec(space_dim) = length(Pop);
Pop_reshape = reshape(Pop,reshape_vec);

repmat_vec = [size_psi(1:(space_dim-1)),1,size_psi((space_dim+1):end)];
Pop_rep = repmat(Pop_reshape,repmat_vec);

psi_TOF=ifft(Pop_rep.*fft(psi,[],space_dim),[],space_dim);

end