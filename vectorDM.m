close all
clear all
clc

%% Philip Mocz (2021), Princeton University
% Merge solitons with Vector Dark Matter Formulation


% Internal units:
% [L] = kpc
% [M] = Msun
% [E] = Msun (km/s)^2

% scaling
% {x, t, rho, m} --> {a x, b t, b^-2 rho, a^-2 b m}

% v = (hbar / m) * grad(phase)


addpath('helpers/')

%% Parameters
m22      = 1;                              % (m/ 10^-22 eV)
Lbox     = 20;                             % kpc
N        = 128;   64;128;256;512; %                       % resolution
Tfinal   = 2;                              % kpc/(km/s) ~ 978 Myr
Nout     = 20; 200;                            % number of output
myseed   = 42; %
saveEnergies = true;     %true;  false;
runScalarVersion = false;  % false; true


output_tag = '';
if runScalarVersion
    output_tag = '_scalar';
end

output   = ['output/vdm_s' num2str(myseed) 'r' num2str(N) 'o' num2str(Nout) output_tag '/'];
energyFile = [output 'energy.txt'];
plotRealTime = false; true; false;


cc = 1;
dt_energySave = 0.001;


%% Constants
hbar = 1.71818131e-87;       % hbar / (mass of sun * (km/s) * kpc)
m = m22 * 8.96215327e-89;    % 10^-22 eV / c^2 / mass of sun
G = 4.3022682e-6;            % G/((km/s)^2*kpc/mass of sun)
c = 299792.458;              % c / (km/s)


m_per_hbar = m/hbar;

clim = [5 9];

%% setup


% spacing & grid
dx = Lbox / N;
xlin = ((0:N-1)' + 0.5) * dx;   % Note, x=0 and x=Lbox are the same point!

%assert(dx < 1); % Make sure we have better than 1 kpc resolution!

[x, y, z] = meshgrid(xlin, xlin, xlin);



% IC
rng(myseed);
Ncores_per_dim = 4;
rc = .4 + .6*rand(3*Ncores_per_dim,1); % \in [.4,1]
ctr = Lbox * rand(3*Ncores_per_dim,3);

psi1 = zeros(size(x));
psi2 = zeros(size(x));
psi3 = zeros(size(x));

for i = 1:3*Ncores_per_dim
    for j = 1:N
        if i <= Ncores_per_dim
            psi1(:,:,j) = psi1(:,:,j) + sqrt(solitonProfile( sqrt( ...
                periodicDist(x(:,:,j),ctr(i,1),Lbox).^2 + ...
                periodicDist(y(:,:,j),ctr(i,2),Lbox).^2 + ...
                periodicDist(z(:,:,j),ctr(i,3),Lbox).^2 ) ...
                , rc(i), m22 ));
        elseif i <= 2*Ncores_per_dim
            psi2(:,:,j) = psi2(:,:,j) + sqrt(solitonProfile( sqrt( ...
                periodicDist(x(:,:,j),ctr(i,1),Lbox).^2 + ...
                periodicDist(y(:,:,j),ctr(i,2),Lbox).^2 + ...
                periodicDist(z(:,:,j),ctr(i,3),Lbox).^2 ) ...
                , rc(i), m22 ));
        else
            psi3(:,:,j) = psi3(:,:,j) + sqrt(solitonProfile( sqrt( ...
                periodicDist(x(:,:,j),ctr(i,1),Lbox).^2 + ...
                periodicDist(y(:,:,j),ctr(i,2),Lbox).^2 + ...
                periodicDist(z(:,:,j),ctr(i,3),Lbox).^2 ) ...
                , rc(i), m22 ));
        end
    end
end

if runScalarVersion
    % in this mode, put all the dark matter into psi1
    psi1 = psi1 + psi2 + psi3;
    psi2 = 0*psi2;
    psi3 = 0*psi3;
end


Mtot = sum(abs(psi1(:)).^2 + abs(psi2(:)).^2 + abs(psi3(:)).^2) * dx^3;       % total mass
rhobar = Mtot / Lbox^3;                  % average density




%stop


clear xlin;
clear x;
clear y;
clear z;



%%


% fourier space variables
fftw('planner','measure');
klin = (-N/2:N/2-1)' * (2*pi/Lbox);
[kx, ky, kz] = meshgrid(klin, klin, klin);
kSq = fftshift(kx.^2 + ky.^2 + kz.^2);
clear klin;
clear kx;
clear ky;
clear kz;



V = zeros(size(psi1));

t = 0;
dt = 0;


% Load any previously saved data, if any
loadedPreviousSave = false;
for snapnum = Nout:-1:0
    filename = [output 'snap' sprintf('%.04d',snapnum) '.h5'];
    if exist(filename,'file')
        t = hdf5read(filename, '/time');
        psi1 = hdf5read(filename, '/psi1Re') + 1.i * hdf5read(filename, '/psi1Im');
        psi2 = hdf5read(filename, '/psi2Re') + 1.i * hdf5read(filename, '/psi2Im');
        psi3 = hdf5read(filename, '/psi3Re') + 1.i * hdf5read(filename, '/psi3Im');
        loadedPreviousSave = true;
        break
    end
end

% else, save first snapshot
if ~loadedPreviousSave
    if ~exist(output,'dir')
        mkdir(output);
    end
    snapnum = 0;
    filename = [output 'snap' sprintf('%.04d',snapnum) '.h5'];
    % if exist(filename,'file')
    %     delete(filename);
    % end
    hdf5write(filename, '/time', double(t))
    hdf5write(filename, '/m22', double(m22), 'WriteMode', 'append')
    hdf5write(filename, '/Lbox', double(Lbox), 'WriteMode', 'append')
    hdf5write(filename, '/psi1Re', double(real(psi1)), 'WriteMode', 'append')
    hdf5write(filename, '/psi1Im', double(imag(psi1)), 'WriteMode', 'append')
    hdf5write(filename, '/psi2Re', double(real(psi2)), 'WriteMode', 'append')
    hdf5write(filename, '/psi2Im', double(imag(psi2)), 'WriteMode', 'append')
    hdf5write(filename, '/psi3Re', double(real(psi3)), 'WriteMode', 'append')
    hdf5write(filename, '/psi3Im', double(imag(psi3)), 'WriteMode', 'append')
    
    % Save Image
    savname = ['frames/s' num2str(myseed) 'r' num2str(N) 'o' num2str(Nout) 's' num2str(snapnum) '.png'];
    A = log10([ mean(abs(psi1).^2,3) mean(abs(psi2).^2,3) mean(abs(psi3).^2,3) ]);
    A = circshift(A, [N/2-round(108*N/256)+1, N/2-round(12*N/256)+1]);
    Amax = clim(2);
    Amin = clim(1);
    A = (A - Amin)/(Amax-Amin);
    A(A<0) = 0;
    A(A>1) = 1;
    Ncolor = 256;
    A = uint8(A*Ncolor);
    mymap = inferno(Ncolor);
    imwrite(A, mymap, savname,'png');
end

snapnum = snapnum + 1;
saveNextTurn = 0;

clear Ncores;
clear rc;
clear ctr;

if plotRealTime
    fh = figure(1);
    set(fh,'position',[0,0,1200,500]);
end

%stop


%%
energyFileID = fopen(energyFile, 'a');


%% simulation
tic;
while t < Tfinal
    
    % potential - (1/2) kick
    psi1 = exp(-1.i * dt/2 * m_per_hbar * V).*psi1;
    psi2 = exp(-1.i * dt/2 * m_per_hbar * V).*psi2;
    psi3 = exp(-1.i * dt/2 * m_per_hbar * V).*psi3;
    
    % kinetic - drift
    psi1 = fftn(psi1);
    psi1 = exp(dt * (1/m_per_hbar)/2 * -1.i * kSq) .*psi1;
    psi1 = ifftn(psi1);
    
    psi2 = fftn(psi2);
    psi2 = exp(dt * (1/m_per_hbar)/2 * -1.i * kSq) .*psi2;
    psi2 = ifftn(psi2);
    
    psi3 = fftn(psi3);
    psi3 = exp(dt * (1/m_per_hbar)/2 * -1.i * kSq) .*psi3;
    psi3 = ifftn(psi3);
    
    % potential - (1/2) kick
    V = 4*pi*G * (abs(psi1).^2 + abs(psi2).^2 + abs(psi3).^2 - rhobar);
    V = -fftn(V);
    V = V ./ ( kSq  + (kSq==0) );
    V = ifftn(V);
    
    psi1 = exp(-1.i * dt/2 * m_per_hbar * V).*psi1;
    psi2 = exp(-1.i * dt/2 * m_per_hbar * V).*psi2;
    psi3 = exp(-1.i * dt/2 * m_per_hbar * V).*psi3;
    
    
    t = t + dt
    
    dt = min( m_per_hbar/6*dx^2, 1./(m_per_hbar*max(abs(V(:)))) );
    
    
    if t+dt >= snapnum * Tfinal/Nout
        dt = snapnum * Tfinal/Nout - t;
        assert(dt > 0);
        saveNextTurn = 1;
    end
    
    % save snapshot
    if saveNextTurn %t > snapnum * Tfinal / Nout
        filename = [output 'snap' sprintf('%.04d',snapnum) '.h5'];
        %         if exist(filename,'file')
        %             delete(filename);
        %         end
        hdf5write(filename, '/time', double(t))
        hdf5write(filename, '/m22', double(m22), 'WriteMode', 'append')
        hdf5write(filename, '/Lbox', double(Lbox), 'WriteMode', 'append')
        hdf5write(filename, '/psi1Re', double(real(psi1)), 'WriteMode', 'append')
        hdf5write(filename, '/psi1Im', double(imag(psi1)), 'WriteMode', 'append')
        hdf5write(filename, '/psi2Re', double(real(psi2)), 'WriteMode', 'append')
        hdf5write(filename, '/psi2Im', double(imag(psi2)), 'WriteMode', 'append')
        hdf5write(filename, '/psi3Re', double(real(psi3)), 'WriteMode', 'append')
        hdf5write(filename, '/psi3Im', double(imag(psi3)), 'WriteMode', 'append')
        
        % Save Image
        savname = ['frames/s' num2str(myseed) 'r' num2str(N) 'o' num2str(Nout) 's' num2str(snapnum) '.png'];
        A = log10([ mean(abs(psi1).^2,3) mean(abs(psi2).^2,3) mean(abs(psi3).^2,3) ]);
        A = circshift(A, [N/2-round(108*N/256)+1, N/2-round(12*N/256)+1]);
        Amax = clim(2);
        Amin = clim(1);
        A = (A - Amin)/(Amax-Amin);
        A(A<0) = 0;
        A(A>1) = 1;
        Ncolor = 256;
        A = uint8(A*Ncolor);
        mymap = inferno(Ncolor);
        imwrite(A, mymap, savname,'png');
        
        
        snapnum = snapnum + 1;
        saveNextTurn = 0;
        
        %% plot
        if plotRealTime
            figure(fh);
            imagesc(log10([  mean(abs(psi1).^2,3) mean(abs(psi2).^2,3)  mean(abs(psi3).^2,3) ]))
            caxis(clim)
            axis off
            pbaspect(gca,[3 1 1])
            colormap(inferno)
            pause(0.0001);
        end
    end
    
    % save enrergies
    if saveEnergies
        if t > dt_energySave * cc
            rho = abs(psi1).^2 + abs(psi2).^2 + abs(psi3).^2;
            V = getpotential( rho, rhobar, Lbox, G );
            [ Krho_1, Kv_1, W_1, KQ_1, ~, ~, ~, ~ ] = getenergies( psi1, V, Lbox, G, m, hbar );
            [ Krho_2, Kv_2, W_2, KQ_2, ~, ~, ~, ~ ] = getenergies( psi2, V, Lbox, G, m, hbar );
            [ Krho_3, Kv_3, W_3, KQ_3, ~, ~, ~, ~ ] = getenergies( psi3, V, Lbox, G, m, hbar );
            fprintf(energyFileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f \n', [t Krho_1, Kv_1, W_1, KQ_1, Krho_2, Kv_2, W_2, KQ_2, Krho_3, Kv_3, W_3, KQ_3 ]);
            cc = cc + 1;
        end
    end
    
    
end
toc;


%%
fclose(energyFileID);


