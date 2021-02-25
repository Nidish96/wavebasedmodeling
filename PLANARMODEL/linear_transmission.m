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

% Lrem = range(INm.X)*80; 
Lrem = 12.5;
% Nerem = fix(Lrem/Le);
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
    maxom_tr2 = sqrt(pars.G/pars.rho)/Lesp*pars.A;
    
%     maxom_tr = maxom_tr*(pars.rho*pars.A);
%     maxom_ax = maxom_ax*(pars.rho*pars.A);
%     maxom_tr2 = maxom_tr2*(pars.rho*pars.Iz);
    
    % BMsp1
    Xqps1 = Nsf*BMsp1.X(n1:n2);
    etax = abs(((Xqps1-Xs1)/Ltots1)).^m;
    
    Cs1(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) = Cs1(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) + ...
        (Nsf'*diag((maxom_ax/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For axial waves
    Cs1(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) = Cs1(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) + ...
        (Nsf'*diag((maxom_tr/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For transverse waves
    Cs1(([n1 n2]-1)*3+3, ([n1 n2]-1)*3+3) = Cs1(([n1 n2]-1)*3+3, ([n1 n2]-1)*3+3) + ...
        (Nsf'*diag((maxom_tr2/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For transverse waves 2
    
%     figure(2)
%     plot(Xqps1, etax, '.-'); hold on
    
    % BMsp2
    Xqps2 = Nsf*BMsp2.X(n1:n2);
    etax = abs(((Xqps2-Xs2)/Ltots2)).^m;
    
    Cs2(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) = Cs2(([n1 n2]-1)*3+1, ([n1 n2]-1)*3+1) + ...
        (Nsf'*diag((maxom_ax/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For axial waves
    Cs2(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) = Cs2(([n1 n2]-1)*3+2, ([n1 n2]-1)*3+2) + ...
        (Nsf'*diag((maxom_tr/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For transverse waves
    Cs2(([n1 n2]-1)*3+3, ([n1 n2]-1)*3+3) = Cs2(([n1 n2]-1)*3+3, ([n1 n2]-1)*3+3) + ...
        (Nsf'*diag((maxom_tr2/Lesp)*etax.*wi)*Nsf)*Lesp/2*pars.A;  % For transverse waves 2
    
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
NsensperBM = 20;  % Number of sensors per BM1/BM2 regions;
SLocs1 = linspace(BM1.X(1), BM1.X(end), NsensperBM+2);  SLocs1([1 end]) = [];
SLocs2 = linspace(BM2.X(1), BM2.X(end), NsensperBM+2);  SLocs2([1 end]) = [];

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
fin = ceil(NsensperBM/2);
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

U0 = zeros(MDL.Ndofs, 1);
Ud0 = zeros(MDL.Ndofs, 1);

opts = struct('Display', 'waitbar');

tic
[T, U, Ud, Udd] = MDL.HHTAMARCH(T0, T1, 1/fsamp, ...
    U0, Ud0, Fex.*(ft(:).*wndw)', opts);
toc

%% Process Data
udd_s1 = Lsens1*Udd;
udd_s2 = Lsens2*Udd;

puls_s1 = zeros(NsensperBM*3, length(T));
puls_s2 = zeros(NsensperBM*3, length(T));

% npr = 3;
% npc = 5;
% figure(3)
% clf()
% 
% figure(4)
% clf()
for si=1:NsensperBM
    for di=1:3
        [s1,fv,t] = STFT([zeros(1, Np2-1), udd_s1((si-1)*3+di,:)], hanning(Np2), 1, floor(length(T))/2, fsamp);        
        puls_s1((si-1)*3+di, :) = interp1(fv, s1, fpuls)/(length(T)/2);
        
        [s2,fv,t] = STFT([zeros(1, Np2-1), udd_s2((si-1)*3+di,:)], hanning(Np2), 1, floor(length(T))/2, fsamp);        
        puls_s2((si-1)*3+di, :) = interp1(fv, s2, fpuls)/(length(T)/2);
        
        
%         figure(3)
%         subplot(npc, npr, (si-1)*3+di)
% %         imagesc(t, fv, abs(s1));
%         imagesc(t, fv, 20*log10(abs(s1)));
% %         xlim([0 0.1])
%         colorbar;
%         set(gca, 'YDir', 'normal')
%         title(sprintf('BM1 Sensor %d: DOF %d', si, di))
%         
%         figure(4)
%         subplot(npc, npr, (si-1)*3+di)
% %         imagesc(t, fv, abs(s2));
%         imagesc(t, fv, 20*log10(abs(s2)));
% %         xlim([0 0.1])
%         colorbar;
%         set(gca, 'YDir', 'normal')
%         title(sprintf('BM2 Sensor %d: DOF %d', si, di))
    end
end

%% Estimate transmission Coefficients
ks = linspace(0, 2*2*pi*32/120e-3, 1000); ks(1) = [];
[ws, wbyks, dwdks] = TIMOWSPDS(ks, pars);  % Solution of dispersion relationship
caxl = sqrt(pars.E/pars.rho);  % Axial wave speed
ws = [ws; caxl*ks];
wbyks = [wbyks; caxl*ones(size(ks))];
dwdks = [dwdks; caxl*ones(size(ks))];

kpuls_t1 = interp1(ws(3,:), ks, 2*pi*fpuls);
kpuls_t2 = interp1(ws(4,:), ks, 2*pi*fpuls);
kpuls_l = interp1(ws(5,:), ks, 2*pi*fpuls);

tmwss = TIMOWSPDS(2*pi*fpuls, pars, 'w');
tmws = [0; 0]*1j;
k = 1;
for i=1:4
    if tmwss(i)>0
        if isreal(tmwss(i))
            tmws(k) = tmwss(i);
            k = k+1;
        elseif imag(tmwss(i))>0
            tmws(k) = tmwss(i);
            k = k+1;
        end
    end
end

[~, m1] = min(abs(tmws-kpuls_t1));
kpuls_t1 = tmws(m1);
m2 = setdiff(1:2, m1);
kpuls_t2 = tmws(m2);

kpuls_l = 2*pi*fpuls/caxl;

% We now have [kpuls_t1, kpuls_t2, kpuls_l]

TMcoefs1_num = zeros(1, NsensperBM*3)*1j;
TMcoefs1_thr = zeros(1, NsensperBM*3)*1j;
[pks, pkin] = findpeaks(abs(puls_s1((fin-1)*3+din, :)));
[~, mi] = max(abs(pks));
pkin = pkin(mi);
for si=1:NsensperBM
    for di=1:3
        % BM1
        [pks, plocs1] = findpeaks(abs(puls_s1((si-1)*3+di, :)));
        [~, mi] = max(abs(pks));
        plocs1 = plocs1(mi);
        
        del_t = T(plocs1) - T(pkin);
        del_x = SLocs1(si) - SLocs1(fin);
        
        TMcoefs1_num(1, (si-1)*3+di) = exp(-1j*del_t*2*pi*fpuls)*puls_s1((si-1)*3+di, plocs1)/puls_s1((fin-1)*3+din, pkin);
        TMcoefs1_thr(1, (si-1)*3+di) = exp(1j*kpuls_t1*del_x);
    end
end

%% All w.r.t SLocs2(1)

TMcoefs2_num = zeros(1, NsensperBM*3)*1j;
TMcoefs2_thr = zeros(1, NsensperBM*3)*1j;
[pks, pkin] = findpeaks(abs(puls_s2((1-1)*3+2, :)));
% [pks, pkin] = findpeaks(abs(puls_s1((fin-1)*3+din, :)));
[~, mi] = max(abs(pks));
pkin = pkin(mi);
for si=1:NsensperBM
    for di=1:3
        % BM2
        [pks, plocs2] = findpeaks(abs(puls_s2((si-1)*3+di, :)));
        [~, mi] = max(abs(pks));
        plocs2 = plocs2(mi);
        
        del_t = T(plocs2) - T(pkin);
        del_x = SLocs2(si) - SLocs2(1);
%         del_x = SLocs2(si) - SLocs1(fin);
        
        TMcoefs2_num(1, (si-1)*3+di) = exp(-1j*del_t*2*pi*fpuls)*puls_s2((si-1)*3+di, plocs2)/puls_s2((1-1)*3+2, pkin);
%         TMcoefs2_num(1, (si-1)*3+di) = exp(-1j*del_t*2*pi*fpuls)*puls_s2((si-1)*3+di, plocs2)/puls_s1((fin-1)*3+din, pkin);
        TMcoefs2_thr(1, (si-1)*3+di) = exp(1j*kpuls_t1*del_x);
    end
end

return
%%
figure(50)
clf()
set(gcf, 'Color', 'white')

% plot(SLocs1, abs(TMcoefs1_num(2:3:end)), 'bo', 'MarkerFaceColor', 'b'); hold on
% plot(SLocs1, abs(TMcoefs1_thr(2:3:end))/2, 'r*')

plot(SLocs2, abs(TMcoefs2_num(2:3:end)), 'bo', 'MarkerFaceColor', 'b'); hold on
plot(SLocs2, abs(TMcoefs2_thr(2:3:end)), 'r*')

legend('Numerical Estimates', 'Theoretical Estimates', 'Location', 'Best')

xlabel('Sensor Locations')
ylabel('Transmission Coefficient (mag)')

export_fig('./FIGS/FPULS_TMsamp.eps', '-depsc')
%%
figure(5)
clf()
set(gcf, 'Color', 'white')

Sd = fix(linspace(0, NsensperBM-1, 4))+1;
Ed = Sd(2:end)-1;  Ed(end) = Ed(end)+1;
Sd = Sd(1:end-1);
NnDd = Ed-Sd+1;
for di=1:3
    subplot(2,3, di)
    aa = plot(T, 20*log10(abs(puls_s2(di:3:end,:))));
    legend(aa(Sd(di):Ed(di)), [repmat('Sensor ', NnDd(di),1) num2str((Sd(di):Ed(di))')], ...
        'Location', 'northeast')
    xlabel('Time (s)')
    ylabel('Acceleration Amplitude (dB)')
    title(sprintf('DOF %d: Red Region', di))
    
    subplot(2,3, 3+di)
	bb = plot(T, 20*log10(abs(puls_s2(di:3:end,:))));
    legend(bb(Sd(di):Ed(di)), [repmat('Sensor ', NnDd(di),1) num2str((Sd(di):Ed(di))')], ...
        'Location', 'northeast')
    xlabel('Time (s)')
    ylabel('Acceleration Amplitude (dB)')
    title(sprintf('DOF %d: Purple Region', di))
    
end
export_fig('./FIGS/FPULSERESP_samppuls.eps', '-depsc')

%% 
figure(10)
clf()
set(gcf, 'Color', 'white')
for i=1:length(MODBMs)
    plot(MODBMs{i}.X, MODBMs{i}.X*0+i, '-', 'LineWidth', 4); hold on
end
plot(SLocs1, SLocs1*0+2, 'kx')
plot(SLocs2, SLocs2*0+4, 'kx')
axis equal
xlabel('X Coordinate')

annotation('textarrow', [0.4 0.5], [0.8 0.6], 'String', 'Central Region')
text(mean(MODBMs{2}.X)*1.5, 0, 'Sensors (x)')
text(mean(MODBMs{4}.X)/1.5, 6, 'Sensors (x)')

% drawbrace([MODBMs{1}.X(end) -4], [MODBMs{1}.X(1) -4], 20, 'LineWidth', 1)


tt = linspace(0, 2*pi);
rr = 2;
plot(rr*cos(tt)+mean(MODBMs{3}.X), rr*sin(tt)+3, 'k--')

plot(MODBMs{1}.X(end)*[1 1], ylim, 'k--')
% plot(MODBMs{2}.X(end)*[1 1], ylim, 'k--')
% plot(MODBMs{3}.X(end)*[1 1], ylim, 'k--')
plot(MODBMs{4}.X(end)*[1 1], ylim, 'k--')
print('./FIGS/FPULS_setupsamp.eps', '-depsc')
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
figure(3)
set(gcf, 'Color', 'white')
clf();
ti = 65;
    for i=1:length(MODBMs)
        for di=1:3
            subplot(3,1, di)            
%             plot(SLocs1, SLocs1*0, 'bo', 'MarkerFaceColor', 'w'); hold on
%             plot(SLocs2, SLocs2*0, 'ro', 'MarkerFaceColor', 'w'); hold on
            
%             plot(SLocs1, Lsens1(di:3:end, :)*U(:, ti), 'bo', 'MarkerFaceColor', 'w'); hold on
%             plot(SLocs2, Lsens2(di:3:end, :)*U(:, ti), 'ro', 'MarkerFaceColor', 'w'); hold on

            plot(MODBMs{i}.X, MODBMs{i}.X*0, '-'); hold on
            plot(MODBMs{i}.X, MODBMs{i}.L(di:3:end,:)*U(:, ti), '.-')
            
            plot(SLocs1(fin), 0, 'k*')
            title(sprintf('Frame %d/%d (%e s): Dof %d', ti, length(T), T(ti), di))
            
            grid on
            ylim((umax(di)+eps)*1.2*[-1 1])
        end
    end
subplot(3,1, 3)
xlabel('X Coordinate (m)')
export_fig('./FIGS/FPULSERESP_samp.eps', '-depsc')

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
            
%             plot(SLocs1, Lsens1(di:3:end, :)*U(:, ti), 'bo', 'MarkerFaceColor', 'w'); hold on
%             plot(SLocs2, Lsens2(di:3:end, :)*U(:, ti), 'ro', 'MarkerFaceColor', 'w'); hold on
            
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