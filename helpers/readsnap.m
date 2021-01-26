function [ t, m22, Lbox, N, psi ] = readsnap( snapdir, snapnum )
%READSNAP read snapshot
%   Detailed explanation goes here

t = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/time');
m22 = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/m22');
%spin = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/spin');
Lbox = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/Lbox');
psi = hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psiRe');
psi = psi + 1.i * hdf5read([snapdir 'snap' sprintf('%.04d',snapnum) '.h5'], '/psiIm');
N = size(psi,1);

end

