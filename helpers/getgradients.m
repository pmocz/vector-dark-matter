function [ dfdx, dfdy, dfdz ] = getgradients( f, Lbox )
%GETGRADIENTS compute gradient of a field
%

N = size(f,1);
dx = Lbox/N;

%finite difference approach

dfdx = (circshift(f,[0 -1 0]) - circshift(f,[0 1 0])) / (2*dx);
dfdy = (circshift(f,[-1 0 0]) - circshift(f,[1 0 0])) / (2*dx);
dfdz = (circshift(f,[0 0 -1]) - circshift(f,[0 0 1])) / (2*dx);


end

