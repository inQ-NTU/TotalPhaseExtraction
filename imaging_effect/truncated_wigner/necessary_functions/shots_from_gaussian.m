function shots = shots_from_gaussian(kernel,no_real)

no_gridpoints = size(kernel,1);

[U,D] = eig(kernel);
eigval = diag(D);
std_modes = sqrt(1./eigval);

coeff = randn(no_real,no_gridpoints);
coeff = coeff*diag(std_modes);

shots = U*(coeff.');
shots = shots';

end