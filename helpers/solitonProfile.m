function rhoc = solitonProfile( r, rc, m22 )
%SOLITONPROFILE returns soliton density at radius r
%   units: kpc, Msun, km/s
%   input: sample radius r (kpc), core radius rc (kpc), m22
%   total mass is ~ 2.20e8 /rc / m22^2

rho0 = 1.9e7 * m22^-2 * rc^-4;    % Msun/kpc^3
rhoc = rho0 ./ (1 + 0.091 * (r/rc).^2).^8;

end

