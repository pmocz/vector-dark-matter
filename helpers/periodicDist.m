function d = periodicDist( x, origin, Lbox )
%PERIODICDIST return minimum distance to origin in periodic domain
%   Detailed explanation goes here

d = x - origin;
d = d + Lbox .* ( d <= Lbox/2 );
d = d - Lbox .* ( d >  Lbox/2 );

end

