close all
clear all
clc
% Philip Mocz (2021), Princeton University
% Do a simple analysis of simulations
% look at energies, resolution

% Internal units:
% [L] = kpc
% [M] = Msun
% [E] = Msun (km/s)^2


%%
Ns = [ 32 64 128 ];

runScalarVersion = true;  % false; true


output_tag = '';
if runScalarVersion
    output_tag = '_scalar';
end

%stop

%% simulation parameters
m22      = 1;                              % (m/ 10^-22 eV)
Lbox     = 20;                             % kpc
%N        = 64; %                          % resolution
Tfinal   = 2;                              % kpc/(km/s) ~ 978 Myr
Nout     = 20;                             % number of output
myseed   = 42;                             % seed


output_root = '../';


addpath('../helpers/')

% constants
hbar = 1.71818131e-87;
G = 4.3022682e-6;


%% Plot Projections
snaps = 0:Nout;
%fh = figure;
%set(fh,'position',[0 0 600 600],'PaperPosition',[0 0 6 6]);

cc = 1;
for N = Ns
    
    
    
    snapdir   = [output_root 'output/vdm_s' num2str(myseed) 'r' num2str(N) 'o' num2str(Nout) output_tag '/'];
    
    
    for snapnum = snaps
        [N snapnum]
        
        
        [ t, m22, Lbox, N, psi1, psi2, psi3 ] = readsnap( snapdir, snapnum );
        m = m22 * 8.96215327e-89;
        
        
        
        
        if snapnum == Nout
            %%  Density Projections
            cmap = inferno(256);
            clim = [5 9];
            % Psi 1
            savname = ['../writeup/s' num2str(myseed) 'r' num2str(N) 's' num2str(snapnum) 'psi1' output_tag '.png'];
            A = log10(mean(abs(psi1).^2,3));
            my_imwrite(A,cmap,clim,savname)
            
            % Psi 2
            savname = ['../writeup/s' num2str(myseed) 'r' num2str(N) 's' num2str(snapnum) 'psi2' output_tag '.png'];
            A = log10(mean(abs(psi2).^2,3));
            my_imwrite(A,cmap,clim,savname)
            
            % Psi 3
            savname = ['../writeup/s' num2str(myseed) 'r' num2str(N) 's' num2str(snapnum) 'psi3' output_tag '.png'];
            A = log10(mean(abs(psi3).^2,3));
            my_imwrite(A,cmap,clim,savname)
            
            
            %% Slice of the Phase
            cmap = redblue(256);
            clim = [-1 1];
            % psi 1
            savname = ['../writeup/s' num2str(myseed) 'r' num2str(N) 's' num2str(snapnum) 'psi1_phase' output_tag '.png'];
            A = cos(angle(psi1(:,:,end/2)));
            my_imwrite(A,cmap,clim,savname)
            
            % psi 2
            savname = ['../writeup/s' num2str(myseed) 'r' num2str(N) 's' num2str(snapnum) 'psi2_phase' output_tag '.png'];
            A = cos(angle(psi2(:,:,end/2)));
            my_imwrite(A,cmap,clim,savname)
            
            % psi 3
            savname = ['../writeup/s' num2str(myseed) 'r' num2str(N) 's' num2str(snapnum) 'psi3_phase' output_tag '.png'];
            A = cos(angle(psi3(:,:,end/2)));
            my_imwrite(A,cmap,clim,savname)
            
        end
        
        
        
        %% Calculate and save energies
        rho = abs(psi1).^2 + abs(psi2).^2 + abs(psi3).^2;
        rhobar = mean( rho(:) );
        V = getpotential( rho, rhobar, Lbox, G );
        [ Krho_1{snapnum+1,N}, Kv_1{snapnum+1,N}, W_1{snapnum+1,N}, KQ_1{snapnum+1,N}, ~, ~, ~, ~ ] = getenergies( psi1, V, Lbox, G, m, hbar );
        [ Krho_2{snapnum+1,N}, Kv_2{snapnum+1,N}, W_2{snapnum+1,N}, KQ_2{snapnum+1,N}, ~, ~, ~, ~ ] = getenergies( psi2, V, Lbox, G, m, hbar );
        [ Krho_3{snapnum+1,N}, Kv_3{snapnum+1,N}, W_3{snapnum+1,N}, KQ_3{snapnum+1,N}, ~, ~, ~, ~ ] = getenergies( psi3, V, Lbox, G, m, hbar );
        
        
        
    end
    
    cc = cc + 1;
    
end




%% Plot energies
fig = figure;

my_colors = lines(5);

tt = snaps * Tfinal/Nout;

for N = fliplr(Ns)
    
    
    for snapnum = snaps
        Krho(snapnum+1) = Krho_1{snapnum+1,N} + Krho_2{snapnum+1,N} + Krho_3{snapnum+1,N};
        Kv(snapnum+1) = Kv_1{snapnum+1,N} + Kv_2{snapnum+1,N} + Kv_3{snapnum+1,N};
        KQ(snapnum+1) = KQ_1{snapnum+1,N} + KQ_2{snapnum+1,N} + KQ_3{snapnum+1,N};
        W(snapnum+1) = W_1{snapnum+1,N} + W_2{snapnum+1,N} + W_3{snapnum+1,N};
    end
    h1 = plot(tt,W,'o-','linewidth',N/Ns(1),'color',my_colors(1,:));
    set(h1, 'markerfacecolor', get(h1, 'color'));
    hold on
    
    h1 = plot(tt,KQ,'o-','linewidth',N/Ns(1),'color',my_colors(2,:));
    set(h1, 'markerfacecolor', get(h1, 'color'));
    
    h1 = plot(tt,Kv,'o-','linewidth',N/Ns(1),'color',my_colors(3,:));
    set(h1, 'markerfacecolor', get(h1, 'color'));
    
    h1 = plot(tt,Krho,'o-','linewidth',N/Ns(1),'color',my_colors(5,:));
    set(h1, 'markerfacecolor', get(h1, 'color'));
    
    h1 = plot(tt,KQ+W,'o-','linewidth',N/Ns(1),'color',my_colors(4,:));
    set(h1, 'markerfacecolor', get(h1, 'color'));
    
end

title('Convergence Plot: $N=32$ (thinnest), $64$, $128$ (thickest)','interpreter','latex')
xlabel('$t$','interpreter','latex','fontsize',14)
ylabel('energy','interpreter','latex','fontsize',14)
axis([0 Tfinal -1e13 1.5e13])

lh = legend('$W$','$K_Q$','$K_v$','$K_\rho$','$W+K_Q$');
set(lh,'location','northwest','interpreter','latex')

saveas(fig,['../writeup/energies' output_tag '.eps'],'epsc2');





%% Plot energies directly from energy file
fig = figure;
set(fig,'position',[0 0 900 600],'PaperPosition',[0 0 9 6]);

for N = fliplr(Ns)
    % read energy file
    snapdir   = [output_root 'output/vdm_s' num2str(myseed) 'r' num2str(N) 'o' num2str(Nout) output_tag '/'];
    energyFile = [snapdir 'energy.txt'];
    fileID = fopen(energyFile);
    A = fscanf(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
    A = reshape(A,13,numel(A)/13);
    fclose(fileID);
    
    t = A(1,:);
    sKrho_1 = A(2,:);
    sKv_1 = A(3,:);
    sW_1 = A(4,:);
    sKQ_1 = A(5,:);
    sKrho_2 = A(6,:);
    sKv_2 = A(7,:);
    sW_2 = A(8,:);
    sKQ_2 = A(9,:);
    sKrho_3 = A(10,:);
    sKv_3 = A(11,:);
    sW_3 = A(12,:);
    sKQ_3 = A(13,:);
    
    sKv = sKv_1 + sKv_2 + sKv_3;
    sKQ = sKQ_1 + sKQ_2 + sKQ_3;
    sKrho = sKrho_1 + sKrho_2 + sKrho_3;
    sW = sW_1 + sW_2 + sW_3;
    
    plot(t,sW,'linewidth',N/Ns(1),'color',my_colors(1,:));
    hold on
    plot(t,sKQ,'linewidth',N/Ns(1),'color',my_colors(2,:));
    plot(t,sKv,'linewidth',N/Ns(1),'color',my_colors(3,:));
    plot(t,sKrho,'linewidth',N/Ns(1),'color',my_colors(5,:));
    plot(t,sW+sKQ,'linewidth',N/Ns(1),'color',my_colors(4,:));
end

title('Convergence Plot: $N=32$ (thinnest), $64$, $128$ (thickest)','interpreter','latex')
xlabel('$t$','interpreter','latex','fontsize',14)
ylabel('energy','interpreter','latex','fontsize',14)
axis([0 Tfinal -1.5e13 1.5e13])

lh = legend('$W$','$K_Q$','$K_v$','$K_\rho$','$W+K_Q$');
set(lh,'location','northwest','interpreter','latex')

saveas(fig,['../writeup/energies2' output_tag '.eps'],'epsc2');
