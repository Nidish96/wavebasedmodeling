clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM/')

% set(0,'defaultTextInterpreter','tex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',16)

Nein = 8;

load(sprintf('./DATS/INTMATS_%dIN.mat', Nein), 'IN1', 'IN2', 'INm', ...
    'Mintf', 'Kintf', 'Mim', 'Kim', 'Fbolt_intf', 'Qrel', 'Trel', 'Xqps', ...
    'pars', 'parsbolt', 'parsint', 'wdt');

%% Building "Infinite Domains"
Le = INm.X(2)-INm.X(1);

ks = linspace(0, 2*2*pi/Le, 2000); ks(1) = [];
[ws, wbyks, dwdks] = TIMOWSPDS(ks, pars);  % Solution of dispersion relationship

Lrem = range(INm.X)*80; 
Nerem = fix(Lrem/Le);
Nerem = 100;
Lerem = Lrem/Nerem;

Nesp = 40;  % Sponge regions
Lesp = Lerem;
Lsp = Lesp*Nesp;

BM1 = struct('X', linspace(INm.X(1)-Lrem, INm.X(1), Nerem+1)', ...
             'Y', zeros(Nerem+1, 1), ...
             'Wi', wdt, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width
BMsp1 = struct('X', linspace(INm.X(1)-Lrem-Lsp, INm.X(1)-Lrem, Nesp+1)', ...
               'Y', zeros(Nesp+1, 1), ...
               'Wi', wdt, ... % in-plane width 
               'Wo', wdt);  % out-of-plane width

BM2 = struct('X', linspace(IN2.X(end), IN2.X(end)+Lrem, Nerem+1)', ...
             'Y', zeros(Nerem+1, 1), ...
             'Wi', wdt, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width
BMsp2 = struct('X', linspace(IN2.X(end)+Lrem, IN2.X(end)+Lrem+Lsp, Nesp+1)', ...
               'Y', zeros(Nesp+1, 1), ...
               'Wi', wdt, ... % in-plane width 
               'Wo', wdt);  % out-of-plane width
           
%% Matrices for monotonic parts
Mr1 = zeros((Nerem+1)*3);
Kr1 = zeros((Nerem+1)*3);
Mr2 = zeros((Nerem+1)*3);
Kr2 = zeros((Nerem+1)*3);

for e=1:Nerem
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    
    % BM1
    [Me, Ke] = PLANTMBMEMATS(BM1.X(n2)-BM1.X(n1), pars);
    Mr1(is:ie, is:ie) = Mr1(is:ie, is:ie) + Me;
    Kr1(is:ie, is:ie) = Kr1(is:ie, is:ie) + Ke;
    
    % BM2
    [Me, Ke] = PLANTMBMEMATS(BM2.X(n2)-BM2.X(n1), pars);
    Mr2(is:ie, is:ie) = Mr2(is:ie, is:ie) + Me;
    Kr2(is:ie, is:ie) = Kr2(is:ie, is:ie) + Ke;
end

%% Matrices for sponge region
Ms1 = zeros((Nesp+1)*3);
Ks1 = zeros((Nesp+1)*3);
Cs1 = zeros((Nesp+1)*3);
Ms2 = zeros((Nesp+1)*3);
Ks2 = zeros((Nesp+1)*3);
Cs2 = zeros((Nesp+1)*3);

m = 4;
fmin = 10;
fmax = 510;
wind = find((min(ws,[],2)<(2*pi*fmin)) & (max(ws,[],2)>(2*pi*fmax)))';
cmax = max(interp1(ws(wind, :), wbyks(wind, :), 2*pi*[fmin fmax]));  % Wave speeds for bending
caxl = sqrt(pars.E/pars.rho);  % Axial wave speed

% figure(2)
% clf()

% figure(3)
% clf()

Nqps = 10;
[xi, wi] = LGWT(Nqps, -1, 1);
Nsf = [1-xi 1+xi]/2;
Xs1 = BMsp1.X(end);  Xe1 = BMsp1.X(1);  Ltots1 = abs(Xe1-Xs1);
Xs2 = BMsp2.X(1);  Xe2 = BMsp2.X(end);  Ltots2 = abs(Xe2-Xs2);
for e=1:Nesp
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    
    % BMsp1
    [Me, Ke] = PLANTMBMEMATS(BMsp1.X(n2)-BMsp1.X(n1), pars);
    Ms1(is:ie, is:ie) = Ms1(is:ie, is:ie) + Me;
    Ks1(is:ie, is:ie) = Ks1(is:ie, is:ie) + Ke;
    
    % BMsp2
    [Me, Ke] = PLANTMBMEMATS(BMsp2.X(n2)-BMsp2.X(n1), pars);
    Ms2(is:ie, is:ie) = Ms2(is:ie, is:ie) + Me;
    Ks2(is:ie, is:ie) = Ks2(is:ie, is:ie) + Ke;
    
    % Integration of Ramped Damping for sponge region
    maxom_tr = cmax/(Lesp*pars.A);
    maxom_ax = caxl/(Lesp*pars.A);
    
    % BMsp1
    Xqps1 = Nsf*BMsp1.X(n1:n2);
    etax = abs(((Xqps1-Xs1)/Ltots1)).^m;
    
    Cs1(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) = Cs1(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) + ...
        (Nsf'*diag((maxom_ax/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For axial waves
    Cs1(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) = Cs1(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) + ...
        (Nsf'*diag((maxom_tr/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For transverse waves
    
%     figure(2)
%     plot(Xqps1, etax, '.-'); hold on
    
    % BMsp2
    Xqps2 = Nsf*BMsp2.X(n1:n2);
    etax = abs(((Xqps2-Xs2)/Ltots2)).^m;
    
    Cs2(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) = Cs2(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) + ...
        (Nsf'*diag((maxom_ax/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For axial waves
    Cs2(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) = Cs2(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) + ...
        (Nsf'*diag((maxom_tr/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For transverse waves
    
%     figure(3)
%     plot(Xqps2, etax, '.-'); hold on
end

%% Assembly
Mmod = blkdiag(Ms1, Mr1, Mim, Mr2, Ms2);
Kmod = blkdiag(Ks1, Kr1, Kim, Kr2, Ks2);
Cmod = blkdiag(Cs1, zeros(size(Kr1)), zeros(size(Kim)), zeros(size(Kr2)), Cs2);

Lass = eye(((Nesp+1+Nerem+1)*2+(Nein+1))*3);

ndcons = cumsum([Nesp+1, Nerem+1, Nein+1, Nerem+1, Nesp+1]);

for i=1:length(ndcons)
    n = ndcons(i);
    Lass((n-1)*3+(1:3), n*3+(1:3)) = eye(3);
end
inds = (ndcons-1)*3+(1:3)';
Lass(:, inds(:)) = [];

Mmod_a = Lass'*Mmod*Lass;
Kmod_a = Lass'*Kmod*Lass;
Cmod_a = Lass'*Cmod*Lass;

%% Null Transform
[Va, Wa] = eig(Kmod_a, Mmod_a);
[Wa, si] = sort(sqrt(diag(abs(Wa))));
Va = Va(:,si);

Lnull = null(Va(:,1:3)'*Mmod_a);
% Lnull = eye(size(Mmod_a));

Mnt = Lnull'*Mmod_a*Lnull;
Knt = Lnull'*Kmod_a*Lnull;
Cnt = Lnull'*Cmod_a*Lnull;

Lassnt = Lass*Lnull;

%% Assembly bookkeeping

MODBMs = {BMsp1, BM1, INm, BM2, BMsp2};
ndstarts = [0, ndcons(1:end-1)]+1;
for i=1:length(ndstarts)
    ns = ndstarts(i);
    n = length(MODBMs{i}.X);
    MODBMs{i}.ndstart = ns;
    MODBMs{i}.L = Lassnt((ns-1)*3+(1:n*3), :);
end

%% "Sensor" Locations
NsensperBM = 5;  % Number of sensors per BM1/BM2 regions;
SLocs1 = linspace(BM1.X(1), BM1.X(end), NsensperBM);
SLocs2 = linspace(BM2.X(1), BM2.X(end), NsensperBM);

% Sensors located analogously on both, and since they're symmetric, there's
% no need to construct interpolants separately
Lsens = zeros(NsensperBM*3, length(BM1.X)*3);
for si=1:NsensperBM
    fe = find((BM1.X(1:end-1)-SLocs1(si)).*(BM1.X(2:end)-SLocs1(si))<=0);
    fe = fe(1);
    
    Nsf = [BM1.X(fe+1)-SLocs1(si) SLocs1(si)-BM1.X(fe)]/(BM1.X(fe+1)-BM1.X(fe));
    Nsf = [Nsf(1) 0 0 Nsf(2) 0 0;
        0 Nsf(1) 0 0 Nsf(2) 0;
        0 0 1/2 0 0 1/2];
    
    n1 = fe;
    n2 = fe+1;
    Lsens((si-1)*3+(1:3), ((n1-1)*3+1):(n2*3)) = Nsf;
end

Lsens1 = [zeros(NsensperBM*3, length(BMsp1.X)*3) Lsens zeros(NsensperBM*3, length([INm.X; BM2.X; BMsp2.X])*3)]*Lassnt;
Lsens2 = [zeros(NsensperBM*3, length([BMsp1.X; BM1.X; INm.X])*3) Lsens zeros(NsensperBM*3, length(BMsp2.X)*3)]*Lassnt;

% Forcing vector at sensor "sin" in BM1 along DOF "din"
fin = 3;
din = 2;
Fex = Lsens1((fin-1)*3+din, :)';

%% GMDOF Model
MDL = MDOFGEN(sparse(Mnt), sparse(Knt), sparse(Cnt), Lassnt);

%% Transient simulation
fsamp = 2^12;
T0 = 0;
T1 = 0.25;
T = (T0:1/fsamp:T1)';
T = T(1:(2^floor(log2(length(T)))));
T1 = T(end);

famp = 1.0;
fpuls = 350;  % pulse frequency (Hz)
ft = sin(2*pi*fpuls*T)*famp;
Np2 = ceil(fsamp/fpuls)*4+2;
wndw = [hanning(Np2); zeros(length(T)-Np2,1)];

% Forcing
% nin = fix(length(MODBMs{2}.X)/2);
% din = 2;
% Fex = MODBMs{2}.L((nin-1)*3+din, :)';

opts = struct('Display', 'waitbar');

tic
[T, U, Ud, Udd] = MDL.HHTAMARCH(T0, T1, 1/fsamp, ...
    zeros(MDL.Ndofs,1), zeros(MDL.Ndofs,1), Fex.*(ft(:).*wndw)', opts);
toc

%% Process Data
udd_s1 = Lsens1*Udd;
udd_s2 = Lsens2*Udd;

puls_s1 = zeros(NsensperBM*3, length(T));
puls_s2 = zeros(NsensperBM*3, length(T));

npr = 3;
npc = 5;
figure(3)
clf()

figure(4)
clf()
for si=1:NsensperBM
    for di=1:3
        [s1,fv,t] = STFT([zeros(1, Np2-1), udd_s1((si-1)*3+di,:)], hanning(Np2), 1, floor(length(T))/2, fsamp);        
        puls_s1((si-1)*3+di, :) = interp1(fv, s1, fpuls)/(length(T)/2);
        
        [s2,fv,t] = STFT([zeros(1, Np2-1), udd_s2((si-1)*3+di,:)], hanning(Np2), 1, floor(length(T))/2, fsamp);        
        puls_s2((si-1)*3+di, :) = interp1(fv, s2, fpuls)/(length(T)/2);
        
        
        figure(3)
        subplot(npc, npr, (si-1)*3+di)
        imagesc(t, fv, abs(s1));
        colorbar;
        set(gca, 'YDir', 'normal')
        title(sprintf('BM1 Sensor %d: DOF %d', si, di))
        
        figure(4)
        subplot(npc, npr, (si-1)*3+di)
        imagesc(t, fv, abs(s2));
        colorbar;
        set(gca, 'YDir', 'normal')
        title(sprintf('BM2 Sensor %d: DOF %d', si, di))
    end
end

%%
figure(5)
clf()

for di=1:3
    subplot(2,3, di)
    plot(T, 20*log10(abs(puls_s1(di:3:end,:))))
    xlabel('Time (s)')
    ylabel('Acceleration Amplitude (dB)')
    
    subplot(2,3, 3+di)
    plot(T, 20*log10(abs(puls_s1(di:3:end,:))))
	if di==3
        legend([repmat('Sensor ', NsensperBM,1) num2str((1:NsensperBM)')])
    end
    xlabel('Time (s)')
    ylabel('Acceleration Amplitude (dB)')
end

%% Plot
figure(1);
clf()
for di=1:3
    subplot(3,1, di)
    plot(T, Lsens1((fin-1)*3+din, :)*U, '.-'); 
    grid on
    
    title(sprintf('Dof %d', di))
    ylabel('Response')
end
xlabel('Time (s)')

pause(1)
figure(2)
clf()

% umax = [max(max(MODBMs{2}.L((nin-1)*3+1, :)*U)), ...
%     max(max(MODBMs{2}.L((nin-1)*3+2, :)*U)), ...
%     max(max(MODBMs{2}.L((nin-1)*3+3, :)*U))];
umax = [max(max(abs(MODBMs{2}.L(1:3:end, :)*U))), ...
    max(max(abs(MODBMs{2}.L(2:3:end, :)*U))), ...
    max(max(abs(MODBMs{2}.L(3:3:end, :)*U)))];
    

for ti=1:length(T)
    clf();
    for i=1:length(MODBMs)
        for di=1:3
            subplot(3,1, di)
            plot(MODBMs{i}.X, MODBMs{i}.X*0, '-'); hold on
            plot(MODBMs{i}.X, MODBMs{i}.L(di:3:end,:)*U(:, ti), '.-')
            
%             plot(SLocs1, SLocs1*0, 'bo', 'MarkerFaceColor', 'w'); hold on
%             plot(SLocs2, SLocs2*0, 'ro', 'MarkerFaceColor', 'w'); hold on
            
            plot(SLocs1, Lsens1(di:3:end, :)*U(:, ti), 'bo', 'MarkerFaceColor', 'w'); hold on
            plot(SLocs2, Lsens2(di:3:end, :)*U(:, ti), 'ro', 'MarkerFaceColor', 'w'); hold on
            
            title(sprintf('Frame %d/%d: Dof %d', ti, length(T), di))
            
            grid on
            ylim((umax(di)+eps)*1.2*[-1 1])
        end
    end
    
    pause(0.001);
end

return 
%% 
[V, Wss] = eig(Kmod_a, Mmod_a);
[Wss, si] = sort(sqrt(abs(diag(Wss))));
V = V(:, si);

sc = 2e-1; ccms = 'rgbym';alph = 0.6;

mi = 5; 
U = V(:, mi);

% sc = 1e-5;
% U = diag(Cmod_a);
% U(1:3:end) = 0;

figure(4); 
clf();

for i=1:length(MODBMs)
    plot(MODBMs{i}.X, MODBMs{i}.X*0, 'k.-'); hold on
    PLANARBMDEPICT(MODBMs{i}.L*U*sc, MODBMs{i}, ccms(i), alph);
end

grid on
title(sprintf('Mode %d: %f Hz', mi, Wss(mi)/2/pi))

axis equal