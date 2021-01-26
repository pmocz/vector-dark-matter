function [ vx, vy, vz ] = getvelocities( psi, Lbox, m, hbar )
%GETVELOCITIES compute the velocity from the wave-function
%   v = nabla S / m
%  psi = sqrt(rho) exp(i S / hbar)

N = size(psi,1);
dx = Lbox/N;

S_per_hbar = angle(psi);
vx = circshift(S_per_hbar,[0 -1 0]) - circshift(S_per_hbar,[0 1 0]);
vy = circshift(S_per_hbar,[-1 0 0]) - circshift(S_per_hbar,[1 0 0]);
vz = circshift(S_per_hbar,[0 0 -1]) - circshift(S_per_hbar,[0 0 1]);
vx(vx > pi) = vx(vx > pi) - 2*pi;
vx(vx <= -pi) = vx(vx <= -pi) + 2*pi;
vy(vy > pi) = vy(vy > pi) - 2*pi;
vy(vy <= -pi) = vy(vy <= -pi) + 2*pi;
vz(vz > pi) = vz(vz > pi) - 2*pi;
vz(vz <= -pi) = vz(vz <= -pi) + 2*pi;
vx = vx / (2*dx) / m * hbar;
vy = vy / (2*dx) / m * hbar;
vz = vz / (2*dx) / m * hbar;

end

