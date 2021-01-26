function [ Krho, Kv, W, KQ, M, Lx, Ly, Lz ] = getenergies( psi, rhobar, Lbox, G, m, hbar )
%GETENERGIES return gradient, kinetic, and potential energies, and mass of
%configuration
%   Detailed explanation goes here

N = size(psi,1);
dx = Lbox/N;

rho = abs(psi).^2;



M = sum(rho(:)) * dx^3;



V = getpotential( rho, rhobar, Lbox, G );
W = sum(rho(:).*V(:)/2) * dx^3;



[ vx, vy, vz ] = getvelocities( psi, Lbox, m, hbar );
Kv = sum(rho(:).*(vx(:).^2 + vy(:).^2 + vz(:).^2)/2) * dx^3;


[ dsrdx, dsrdy, dsrdz ] = getgradients( sqrt(rho), Lbox );
Krho = hbar^2/m^2/2 * sum(dsrdx(:).^2 + dsrdy(:).^2 + dsrdz(:).^2) * dx^3;

[ dpsidx, dpsidy, dpsidz ] = getgradients( psi, Lbox );
KQ = hbar^2/m^2/2 * sum( abs(dpsidx(:)).^2 + abs(dpsidy(:)).^2 + abs(dpsidz(:)).^2 ) * dx^3;


% ang momentum
    dx = Lbox / N;
    xlin = ((0:N-1)' + 0.5) * dx - 0.5*Lbox;
    [x, y, z] = meshgrid(xlin, xlin, xlin);
    Lx = sum( rho(:) .* ( vz(:).*y(:) - vy(:).*z(:)) ) * dx^3;
    Ly = sum( rho(:) .* ( vx(:).*z(:) - vz(:).*x(:)) ) * dx^3;
    Lz = sum( rho(:) .* ( vy(:).*x(:) - vx(:).*y(:)) ) * dx^3;
    
    clear x;
    clear y;
    clear z;

end

