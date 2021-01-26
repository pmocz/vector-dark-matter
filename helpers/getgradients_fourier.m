function [ dfdx, dfdy, dfdz ] = getgradients_fourier( f, Lbox, real_flag )
%GETGRADIENTS compute gradient of a field
% data must be periodic!!
N = size(f,1);
%dx = Lbox/N;

%finite difference approach

% dfdx = (circshift(f,[0 -1 0]) - circshift(f,[0 1 0])) / (2*dx);
% dfdy = (circshift(f,[-1 0 0]) - circshift(f,[1 0 0])) / (2*dx);
% dfdz = (circshift(f,[0 0 -1]) - circshift(f,[0 0 1])) / (2*dx);

% fourier approach

klin = (2*pi/Lbox) * (-N/2:N/2-1)';
[kx, ky, kz] = meshgrid(klin, klin, klin);
kx = fftshift(kx);
ky = fftshift(ky);
kz = fftshift(kz);
f = fftn(f);
dfdx = ifftn(1.i * kx .* f);
dfdy = ifftn(1.i * ky .* f);
dfdz = ifftn(1.i * kz .* f);

if real_flag
    dfdx = real(dfdx);
    dfdy = real(dfdx);
    dfdz = real(dfdx);
end

end

