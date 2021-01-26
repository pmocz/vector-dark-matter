function [ t, m22, Lbox, N, psi1, psi2, psi3 ] = readsnap( snapdir, snapnum )
%READSNAP read snapshot
%   Detailed explanation goes here

t = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/time');
m22 = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/m22');
Lbox = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/Lbox');
psi1 = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psi1Re');
psi1 = psi1 + 1.i * hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psi1Im');
psi2 = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psi2Re');
psi2 = psi2 + 1.i * hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psi2Im');
psi3 = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psi3Re');
psi3 = psi3 + 1.i * hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psi3Im');
N = size(psi1,1);

end

