function J = redblue( m )
%MYJET colorbar
% (|data|)^alpha*sign(data)

% From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
m1 = m*0.5;
r = (0:m1-1)'/max(m1-1,1);
g = r;
r = [r; ones(m1,1)];
g = [g; flipud(g)];
b = flipud(r);

J = [r g b];


end

