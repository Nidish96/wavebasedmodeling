clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM/')

% set(0,'defaultTextInterpreter','tex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',16)

%% 
wdt = 25.4e-3;
nu = 0.3;
kappa = 10*(1+nu)/(12+11*nu);
pars = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu))*kappa, ...
    'A', wdt^2, 'Iz', wdt^4/12);  % Monolithic region 
L = 10;

Cl = sqrt(pars.E/pars.rho);
Cs = sqrt(pars.G/pars.rho);

Wc = sqrt(pars.G*pars.A/(pars.rho*pars.Iz));
fc = Wc/2/pi;

%% Finite Element
Ne = 50;
Le = L/Ne;
Nesp = 20;
Lsp = Le*Nesp;
Xn = ((-Lsp):(Le):(L+Lsp))';
Netot = Ne+Nesp*2;

[Me, Ke] = PLANTMBMEMATS(Le, pars);

M = zeros((Netot+1)*3);
K = zeros((Netot+1)*3);
C = zeros((Netot+1)*3);
for e=1:Netot
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    
    M(is:ie, is:ie) = M(is:ie, is:ie) + Me;
    K(is:ie, is:ie) = K(is:ie, is:ie) + Ke;
end
% C(1,1)         = Cl*(pars.rho*pars.A);
% C(2,2)         = Cs*(pars.rho*pars.A);
% C(3,3)         = Cl*(pars.rho*pars.Iz);
% 
% C(end-2,end-2) = Cl*(pars.rho*pars.A);
% C(end-1,end-1) = Cs*(pars.rho*pars.A);
% C(end,end)     = Cl*(pars.rho*pars.Iz);
%% Sponge region damping
m = 4;

Nqps = 10;
[xi, wi] = LGWT(Nqps, -1, 1);
xi = xi(end:-1:1);
wi = wi(end:-1:1);
Nsf = [1-xi 1+xi]/2;
Xstart1 = Xn(Nesp+1);
Xstart2 = Xn(Nesp+Ne+1);
for e=1:Nesp
    maxom_l = (pars.rho*pars.A)*Cl;
    maxom_1 = (pars.rho*pars.A)*Cs;
    maxom_2 = (pars.rho*pars.Iz)*Cl;
    
    % Left Half
    n1 = e;
    n2 = e+1;
    
    ix = ([n1 n2]-1)*3+1;
    iy = ([n1 n2]-1)*3+2;
    itz = ([n1 n2]-1)*3+3;
    
    Xs = Nsf*Xn(n1:n2);
    etr_x = maxom_l*(abs(Xs-Xstart1)/Lsp).^m;
    etr_y = maxom_1*(abs(Xs-Xstart1)/Lsp).^m;
    etr_tz = maxom_2*(abs(Xs-Xstart1)/Lsp).^m;
    
    C(ix, ix) = C(ix, ix) + Nsf'*diag(etr_x.*wi*Le/2)*Nsf;
    C(iy, iy) = C(iy, iy) + Nsf'*diag(etr_y.*wi*Le/2)*Nsf;
    C(itz, itz) = C(itz, itz) + Nsf'*diag(etr_tz.*wi*Le/2)*Nsf;
    
%     et1 = etr_x;
%     mx1 = Nsf'*diag(etr_x.*wi*Le/2)*Nsf;
    
    % Right Half
    n1 = (Nesp+Ne)+e;
    n2 = (Nesp+Ne)+e+1;
    
    Xs = Nsf*Xn(n1:n2);
    etr_x = maxom_l*(abs(Xs-Xstart2)/Lsp).^m;
    etr_y = maxom_1*(abs(Xs-Xstart2)/Lsp).^m;
    etr_tz = maxom_2*(abs(Xs-Xstart2)/Lsp).^m;
    
    ix = ([n1 n2]-1)*3+1;
    iy = ([n1 n2]-1)*3+2;
    itz = ([n1 n2]-1)*3+3;
    
    C(ix, ix) = C(ix, ix) + Nsf'*diag(etr_x.*wi*Le/2)*Nsf;
    C(iy, iy) = C(iy, iy) + Nsf'*diag(etr_y.*wi*Le/2)*Nsf;
    C(itz, itz) = C(itz, itz) + Nsf'*diag(etr_tz.*wi*Le/2)*Nsf;
    
%     et2 = etr_x;
%     mx2 = Nsf'*diag(etr_x.*wi*Le/2)*Nsf;
end
% C = C*0;
C(1:3:end, 1:3:end) = 0;
C(1,1)         = Cl*(pars.rho*pars.A);
C(end-2, end-2)= Cl*(pars.rho*pars.A);
%% Mid-Point Recovery Matrix
RECOV = zeros(1, Netot+1);
if mod(Netot, 2)==0
    RECOV(Netot/2+1) = 1;
else
    RECOV([(Netot+1)/2 (Netot+3)/2]) = 1/2;
end

%% Null-Space Transformation
[V, D] = eig(K, M);
[D, si] = sort(diag(D));
V = V(:, si);

Lnull = null(V(:, 1:3)'*M);

%% Calculate Wavenumbers for given frequency
% W = 2*pi*1000;  % Frequency of wave of interest  % for X
W = 2*pi*100;
kl = W/Cl;

ktran = TIMOWSPDS(W, pars, 'w');
k1 = abs(real(ktran(2)));
k2 = ktran(3);

if ~isreal(k1) || k1<0 || isreal(k2) || imag(k2)>0 || abs(abs(angle(k2))-pi/2)>1e-10
    error('Please choose wavenumbers manually')
end
k2 = -imag(k2);

wav.P = k1*(1-(W/k1*Cs)^2);
wav.N = k2*(1+(W/k2*Cs)^2);

%% Setup initial conditions
al = [0.0 0.0];  % a_\ell^+, a_\ell^+
a1 = [1.0 0.0]*1e-3;  % a_1^+, a_1^+
a2 = [0.0 0.0];  % a_2^+, a_2^+

wndw = [zeros(Nesp,1); hanning(Ne+1); zeros(Nesp,1)];

x = (Xn((Nesp+1):(Nesp+Ne+1))-mean(Xn((Nesp+1):(Nesp+Ne+1))));
xL = Xn(Nesp+Ne+3);
wndw = [zeros(Nesp,1); cos(pi*x/xL).^2; zeros(Nesp, 1)];
wndwdot = [zeros(Nesp,1); -sin(2*pi*x/xL)*pi/xL; zeros(Nesp, 1)]*0;

ux0 = [exp(-1i*kl*Xn) exp(1i*kl*Xn)]*al';
uxd0 = (1j*W)*ux0;

uy0 = [exp(-1i*k1*Xn) exp(1i*k1*Xn) exp(-k2*Xn) exp(k2*Xn)]*[a1'; a2'];
uyd0 = (1j*W)*uy0;

thz0 = [-(1i*wav.P)*exp(-1i*k1*Xn) (1i*wav.P)*exp(1i*k1*Xn) -(wav.N)*exp(-k2*Xn) (wav.N)*exp(k2*Xn)]*[a1'; a2'];
thzd0 = (1j*W)*thz0;

ux0 = ux0.*wndw;  uxd0 = uxd0.*wndw + ux0.*wndwdot;
uy0 = uy0.*wndw;  uyd0 = uyd0.*wndw + uy0.*wndwdot;
thz0 = thz0.*wndw;  thzd0 = thzd0.*wndw + thz0.*wndwdot;

U0 = real(reshape([ux0 uy0 thz0]', [], 1));
Ud0 = real(reshape([uxd0 uyd0 thzd0]', [], 1));

%% MDOFGEN MODEL
MDL = MDOFGEN(Lnull'*M*Lnull, Lnull'*K*Lnull, Lnull'*C*Lnull, Lnull);

%% Time Transient Simulation

T0 = 0;
T1 = 0.1;
fsamp = 2^14;  % 2^14, 2^18
dt = 1/fsamp;
T = (T0:dt:T1)';
T = T(1:2^fix(log2(length(T))));
T1 = T(end);
Nt = length(T);

%% Lsim
I = eye(MDL.Ndofs);
O = zeros(MDL.Ndofs);
sys = ss([O I;-MDL.M\MDL.K -MDL.M\MDL.C], zeros(MDL.Ndofs*2, 1), [I O], 0);

tic
[Yl, Tl, Xl] = lsim(sys, zeros(Nt,1), T, [Lnull'*U0; Lnull'*Ud0]);
toc
Yld = Xl(:, size(Yl,2)+1:end);

%% HHTA
opts = struct('Display', 'waitbar');
tic
[T, U, Ud, Udd] = MDL.HHTAMARCH(T0, T1, dt, Lnull'*U0, Lnull'*Ud0, zeros(MDL.Ndofs, Nt), opts);
toc
Nt = length(T);
%% Analyze mid-point
uxt = RECOV*Lnull(1:3:end, :)*Ud;
uyt = RECOV*Lnull(2:3:end, :)*Ud;
thzt = RECOV*Lnull(3:3:end, :)*Ud;

uxt_l = RECOV*Lnull(1:3:end, :)*Yld';
uyt_l = RECOV*Lnull(2:3:end, :)*Yld';
thzt_l = RECOV*Lnull(3:3:end, :)*Yld';

hwdw = hanning(Nt);
[fs, UXf] = FFTFUN(T, uxt(:).*hwdw);
[fs, UYf] = FFTFUN(T, uyt(:).*hwdw);
[fs, THZf] = FFTFUN(T, thzt(:).*hwdw);

figure(10)
clf()

plot(T, uyt, '-'); grid on; hold on
plot(Tl, uyt_l, '-')
ylim(max(abs(uyt))*[-1 1])

xlabel('Time (s)')
ylabel('Response')
%% Plot
umax = [max(max(abs(Lnull(1:3:end, :)*U))), max(max(abs(Lnull(2:3:end, :)*U))), max(max(abs(Lnull(3:3:end, :)*U)))]+eps;

figure(2)
clf()

for ti=1:Nt
    clf()
    for i=1:3
        subplot(3,1, i)
        % plot(Xn, Lnull(i:3:end, :)*U(:, ti), '.-'); hold on
        plot(Xn, Lnull(i:3:end, :)*Yl(ti, :)', '.-'); hold on
        plot(Xn(1:Nesp+1), Xn(1:Nesp+1)*0, 'k.')
        plot(Xn((Nesp+Ne)+(1:Nesp+1)), Xn(1:Nesp+1)*0, 'k.')
        grid on
        
        ylim(umax(i)*[-1 1]);
        
        ylabel(sprintf('DOF %d', i))
        
        if (i==1)
            title(sprintf('Frame %d/%d', ti, Nt))
        end
    end
	xlabel('X Coordinate')
    
    pause(0.1);
end
