function V = getpotential( rho, rhobar, Lbox, G )
%GETPOTENTIAL return the potential (periodic) of density distribution
%   Detailed explanation goes here

N = size(rho,1);
%dx = Lbox/N;
%rhobar = sum(rho(:)) * dx^3 / Lbox^3;

% fourier space variables
klin = (2*pi/Lbox) * (-N/2:N/2-1)';
[kx, ky, kz] = meshgrid(klin, klin, klin);
kSq = kx.^2 + ky.^2 + kz.^2;
kSq = fftshift(kSq);
clear kx;
clear ky;
clear kz;

% Poisson solver
V = -(fftn( 4*pi*G * (rho-rhobar) )) ./ ( kSq  + (kSq==0));  % (Vhat)
V = ifftn(V);

% normalize so mean potential is 0
V = V - mean(V(:));

end

