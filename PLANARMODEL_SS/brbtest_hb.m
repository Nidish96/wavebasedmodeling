%clc
clear all 
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/SOLVERS')
addpath('../ROUTINES/HARMONIC')
addpath('../ROUTINES/PLANARBM/')

set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxesFontSize',16)

load('./DATS/BRBMATS.mat', 'Mbrb', 'Kbrb', 'Kbolt', 'Mbolt', 'Fbolt', 'Fboltd', 'Qrel', ...
            'Trel', 'BM1', 'IN1', 'BM2', 'IN2', 'Nein', 'Nemono', 'pars', ...
            'parsint', 'Ri1', 'Ri2', 'L1', 'L2', 'wdt', 'Nqp', 'Q1', 'T1', ...
            'Q2', 'T2', 'BoltLocs', 'Lint', 'Lmono', 'nu')
Prestress = 12e3;
Fbolto = Fbolt;
Fbolt = Fboltd;

%% Getting rid of Rigid Body modes 
Lf = null(Ri1([1:3:end 2:3:end], :)*L1 - Ri2([1:3:end 2:3:end], :)*L2);
[Vf, Wf] = eig(Lf'*(Kbrb+Kbolt)*Lf, Lf'*(Mbrb+Mbolt)*Lf);
[Wf, si] = sort(sqrt(diag(Wf)));
Vf = Vf(:, si);

Ln = Lf*null(Vf(:, 1:3)'*(Lf'*(Mbrb+Mbolt)*Lf));
Ln = null((Lf*Vf(:, 1:3))'*(Mbrb+Mbolt));

%% Quadrature Matrices for relative displacement and tractions 
Nqp = 3;
Q1 = zeros((Nein*Nqp)*2, (Nein+1)*3);
Q2 = zeros((Nein*Nqp)*2, (Nein+1)*3);
T1 = zeros((Nein*Nqp)*2, (Nein+1)*3);
T2 = zeros((Nein*Nqp)*2, (Nein+1)*3);
for e=1:Nein
  n1 = e;
  n2 = e+1;
  
  is = (n1-1)*3+1;
  ie = n2*3;
  
  Le = IN1.X(n2)-IN1.X(n1);
  
  [xi, wi] = LGWT(Nqp, 0, Le); 
  wi = wi*IN1.Wo;  % Multiplying with out of plane width 
  for qi=1:Nqp
    Q1((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 -wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2];
    Q2((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2];
    
    T1((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 -wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2]*wi(qi);
    T2((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2]*wi(qi);
  end 
end
Qrel = Q1*L1((Nemono+1)*3+(1:(Nein+1)*3), :) - Q2*L2(1:(Nein+1)*3, :);
Trel = (T1*L1((Nemono+1)*3+(1:(Nein+1)*3), :) - T2*L2(1:(Nein+1)*3, :))';

LTrel = Ln'*Trel;
QrelL = Qrel*Ln;

%% Contact Model
Aint = Lint*wdt;

gap  = 0;
mu   = 0.25;

Pint = Prestress/Aint;
sint = 1e-5;
chi  = 2;
ktkn = chi*(1-nu)/(2-nu);
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn   = kt/ktkn

K0 = Trel*diag(reshape(kron([kt kn], ones(Nqp*Nein,1))', [], 1))*Qrel;

%% Linear Dissipation
[V0, W0] = eig(Ln'*(Kbrb+Kbolt+K0)*Ln, Ln'*(Mbrb+Mbolt)*Ln);
[W0, si] = sort(sqrt(diag(W0)));
V0 = V0(:, si);
V0 = V0./sqrt(diag(V0'*(Ln'*Mbrb*Ln)*V0)');

zetas = [0.1; 0.2]*1e-2;
ab = [1./(2*W0(1:length(zetas))), W0(1:length(zetas))/2]\zetas;

Cbrb = ab(1)*Mbrb + ab(2)*(Kbrb+Kbolt+K0);

%% Input Force Vector 
Finp = zeros((Nemono+Nein+2)*2*3,1);
% Finp(1) = 1;
Finp(2) = 1;

Finp = [L1;L2]'*Finp;

%% GMDOF
MDL = MDOFGEN(Ln'*(Mbrb+Mbolt)*Ln, Ln'*(Kbrb+Kbolt)*Ln, Ln'*Cbrb*Ln, Ln);

fnl = @(t, u, varargin) ELDRYFRICT2D(t, u, kt, kn, mu, gap, varargin{:});
MDL = MDL.SETNLFUN(2+5, QrelL, fnl, LTrel);

%% Static Prestress Analysis
U0 = (MDL.K + Ln'*(K0)*Ln)\(Ln'*Fbolt*Prestress);

opts = struct('reletol', 1e-6, 'Display', true);
[Ustat, ~, ~, ~, J0] = NSOLVE(@(U) MDL.STATRESFUN(U, Ln'*Fbolt*Prestress), U0, opts);
U = Ln*Ustat;

sc = 1e3; ccm1 = 'b';ccm2 = 'r';alph = 0.6;

figure(1)
clf()
PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph);
colorbar('southoutside')
xlim([-0.01 0.73])
ylim([-1 1]*0.125)

%% Prestressed Linear Modes
[Vp, Wp] = eig(J0, MDL.M);
[Wp, si] = sort(sqrt(diag(Wp)));
Vp       = Vp(:, si);

%% HBM - Force Controlled
h = [0 1 2 3];
Nhc = sum(h==0)+2*sum(h~=0);

Nt = 2^8;
Nd = MDL.Ndofs;

fa = 25.00;  % 0.01, 0.20, 1.00, 5.00, 15.00, 25.00, 

Fl = zeros(Nd*Nhc, 1);
Fl(1:Nd) = Ln'*Fbolt*Prestress;
Fl(Nd+(1:Nd)) = fa*Ln'*Finp;

Wst = 250*2*pi;
Wen = 270*2*pi;

% Wst = 250*2*pi;
% Wen = 255*2*pi;

dw = 1.0;
dsmax = 5;
dsmin = 0.1;

E = HARMONICSTIFFNESS(MDL.M, MDL.C, J0, Wst, h);
U0 = E\Fl;
U0(1:Nd) = Ustat;

Copt = struct('Nmax', 400, 'itDisplay', false, 'angopt', 1e-4, ...
    'dsmin', dsmin, 'arclengthparm', 'orthogonal', 'dsmax', dsmax, ...
    'lsrch', 0, 'DynDscale', 0);
% Dscale = ones(size(U0));
% Dscale(1:Nd) = ones(Nd,1)*max(abs(U0(1:Nd)));
% Dscale(Nd+1:end) = ones(Nd*(Nhc-1),1)*max(abs(U0(Nd+1:end)));

% Dscale = ones(size(U0))*max(abs(U0(1:Nd)));
% Dscale = U0;

% Copt.Dscale = [abs(U0); Wst];

UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
% UC = CONTINUE(@(Uw) HBMRESFUN_TMP(Uw, MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

% Umi = NSOLVE(@(U) HBMRESFUN_TMP([U; Um(end)], MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, 1e-6), Um(1:end-1), nopts);

uout_h = (kron(eye(Nhc), Finp'*Ln)*UC(1:end-1,:));
frf = (uout_h(2,:)-uout_h(3,:)*1j)/fa;

figure(10)
% clf(); 

subplot(2,1, 1)
plot(UC(end,:)/2/pi, 20*log10(abs(frf)), '.-'); hold on 
% plot(UC(end,:)/2/pi, 20*log10(sqrt([1*0 0.5*ones(1, Nhc-1)]*(kron(eye(Nhc), Finp'*Ln)*UC(1:end-1,:)).^2)/fa), '.-'); hold on 
grid on
xlabel('Frequency (Hz)')
ylabel('Dynamic Response Amplitude (dB)')

subplot(2,1, 2)
plot(UC(end,:)/2/pi, rad2deg(angle(frf)), '.-'); hold on 
grid on
xlabel('Frequency (Hz)')
ylabel('Response Phase (degs)')

%% HBM - Amplitude Controlled
h = [0 1];
Nhc = sum(h==0)+2*sum(h~=0);

Nt = 2^7;
Nd = MDL.Ndofs;

Amp = 1e-2;  % 1e-6 1e-5 1e-4

Fl = zeros(Nd*Nhc, 1);
Fl(1:Nd) = Ln'*Fbolt*Prestress;
Fl(Nd+(1:Nd)) = Ln'*Finp;

Wst = 245*2*pi;
Wen = 265*2*pi;

% Wst = 255*2*pi;
% Wen = 252*2*pi;

dw = 1.0;
dsmax = 50;
dsmin = 0.1;

E = HARMONICSTIFFNESS(MDL.M, MDL.C, J0, Wst, h);
U0 = E\Fl;
U0(1:Nd) = Ustat;

Copt = struct('Nmax', 1000, 'itDisplay', false, 'angopt', 1e-1, ...
    'dsmin', dsmin, 'arclengthparm', 'orthogonal', 'dsmax', dsmax, ...
    'lsrch', 0, 'DynDscale', 0);
% Dscale = ones(size(U0));
% Dscale(1:Nd) = ones(Nd,1)*max(abs(U0(1:Nd)));
% Dscale(Nd+1:end) = ones(Nd*(Nhc-1),1)*max(abs(U0(Nd+1:end)));

% Dscale = ones(size(U0))*max(abs(U0(1:Nd)));
% Dscale = U0;

% Copt.Dscale = [Dscale; 1.0];

cons = struct('order', 1, 'type', 'D', 'mask', (Ln'*Finp)', 'sc', 1e8);
% MDL.ACHBRESFUN([U0; 0.01; Wst], Amp, Fl, cons, h, Nt, 1e-6)

% Copt.Dscale = [abs(U0); 1; Wst];
UC = CONTINUE(@(Ufw) MDL.ACHBRESFUN(Ufw, Amp, Fl, cons, h, Nt, 1e-6), [U0; 0.01], Wst, Wen, dw, Copt);

uout_h = (kron(eye(Nhc), Finp'*Ln)*UC(1:end-2,:));
frf = (uout_h(2,:)-uout_h(3,:)*1j)./UC(end-1,:);

figure(10)
% clf(); 

subplot(2,1, 1)
plot(UC(end,:)/2/pi, 20*log10(abs(frf)), '+-'); hold on 
% plot(UC(end,:)/2/pi, 20*log10(sqrt([1*0 0.5*ones(1, Nhc-1)]*(kron(eye(Nhc), Finp'*Ln)*UC(1:end-1,:)).^2)/fa), '.-'); hold on 
grid on
xlabel('Frequency (Hz)')
ylabel('Dynamic Response Amplitude (dB)')

subplot(2,1, 2)
plot(UC(end,:)/2/pi, rad2deg(angle(frf)), '.-'); hold on 
grid on
xlabel('Frequency (Hz)')
ylabel('Response Phase (degs)')

%% Step from Last point
WW = UC(end, end);
UG = UC(1:end-1, end);
[R, dRdUf, dRdw] = MDL.ACHBRESFUN([UG; WW], Amp, Fl, cons, h, Nt, 1e-6);

tgt = [-dRdUf\R; 1];
tgt = tgt/norm(tgt);

WW = WW+1;

nopts = struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', 1e-6, 'utol', 1e-6, ...
    'Display', true, 'lsrch', 0)
NSOLVE(@(Uf) MDL.ACHBRESFUN([Uf; WW], Amp, Fl, cons, h, Nt, 1e-6), UG, nopts);

%%
fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter', 'UseParallel', true);
cons = struct('order', 1, 'type', 'D', 'mask', (Ln'*Finp)', 'sc', 1e8);

fsolve(@(Uf) MDL.ACHBRESFUN([Uf; WW+0.1], Amp, Fl, cons, h, Nt, 1e-6), UG, fopts);

%% Check 
U0 = UG;
W0 = WW+0.1;

[R0, dR0] = MDL.ACHBRESFUN([U0; W0], Amp, Fl, cons, h, Nt, 1e-6);
dRnum = zeros(size(dR0));

hm = 1e-6;
parfor hi=1:Nd*Nhc+1
    hv = zeros(Nd*Nhc+1, 1);
    hv(hi) = 1;
    Rp = MDL.ACHBRESFUN([U0+hv*hm;W0], Amp, Fl, cons, h, Nt, 1e-6) ;
    Rm = MDL.ACHBRESFUN([U0-hv*hm;W0], Amp, Fl, cons, h, Nt, 1e-6) ;
    dRnum(:, hi) = (Rp-Rm)/(2*hm);
    hv(hi) = 0;
    
    fprintf('%d/%d\n', hi, Nd*Nhc);
end

hi = 1;

%% NLVIB
Wst = 1655;
Wen = 1620; 

Sopt = struct('jac', 'full', 'dynamicDscale', 1);
UnC = solve_and_continue(U0, @(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), Wst, Wen, dw, Sopt);

%%
Sopt = struct('jac', 'full', 'dynamicDscale', 1);
UnC = solve_and_continue(U0, @(Uw) HBMRESFUN_TMP(Uw, MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, tol), ...
    Wst, Wen, dw, Copt);

%% Check Gradients
rng(1)
U0 = rand(Nd*Nhc, 1);
W0 = 270*2*pi;

Fl = zeros(Nd*Nhc, 1);
Fl(1:Nd) = Ln'*Fbolt*Prestress;
Fl(Nd+(1:Nd)) = fa*Ln'*Finp;

U0 = HARMONICSTIFFNESS(MDL.M, MDL.C, J0, W0, h)\Fl;

[R0, dR0] = MDL.HBRESFUN([U0; W0], Fl, h, Nt, 1e-6);
dRnum = zeros(size(dR0));

hm = 1e-6;
parfor hi=1:Nd*Nhc
    hv = zeros(Nd*Nhc, 1);
    hv(hi) = 1;
    Rp = MDL.HBRESFUN([U0+hv*hm;W0], Fl, h, Nt, 1e-6) ;
    Rm = MDL.HBRESFUN([U0-hv*hm;W0], Fl, h, Nt, 1e-6) ;
    dRnum(:, hi) = (Rp-Rm)/(2*hm);
    hv(hi) = 0;
    
    fprintf('%d/%d\n', hi, Nd*Nhc);
end

disp([reshape(dR0(:, hi),3,[])' reshape(dR0(:,hi)-dRnum(:, hi),3,[])'])

%%
rng(1)
U0 = rand(Nd*Nhc, 1);
W0 = 270*2*pi;

opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter', 'UseParallel', true);
fsolve(@(U) MDL.HBRESFUN([U; W0], Fl, h, Nt, 1e-6), U0, opts);

%%
nopts = struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', 1e-6, 'utol', 1e-6, ...
    'Display', true, 'lsrch', 0)
Usol = NSOLVE(@(U) MDL.HBRESFUN([U; W0], Fl, h, Nt, 1e-6), U0, nopts);

%%
nopts.lsrch = 1;
NSOLVE(@(U) MDL.HBRESFUN([U; 290*2*pi], Fl, h, Nt, 1e-6), Usol, nopts);

%%
fsolve(@(U) MDL.HBRESFUN([U; 290*2*pi], Fl, h, Nt, 1e-6), Usol, opts);

[R, dR] = MDL.HBRESFUN([Usol; 290*2*pi], Fl, h, Nt, 1e-6);

%% 
fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter', 'UseParallel', true);
cons = struct('order', 1, 'type', 'D', 'mask', (Ln'*Finp)');
fsolve(@(Uf) MDL.ACHBRESFUN([Uf; Wst], Amp, Fl, cons, h, Nt, 1e-6), [U0; 0.01], fopts);

%%
nopts = struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', 1e-6, 'utol', 1e-6, ...
    'Display', true, 'lsrch', 0)
NSOLVE(@(Uf) MDL.ACHBRESFUN([Uf; Wst], Amp, Fl, cons, h, Nt, 1e-6), [U0; 0.01], nopts);

%%
rng(1);
U0 = rand(Nd*Nhc, 1);
W0 = 270*2*pi;

tol = 1e-6;

[R, dRU, dRw] = HBMRESFUN_TMP([U0; W0], MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, tol);

fopts = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'iter', 'UseParallel', true);

U0 = HARMONICSTIFFNESS(MDL.M, MDL.C, MDL.K, W0, h)\Fl;
fsolve(@(U) HBMRESFUN_TMP([U; W0], MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, tol), U0, fopts);

NSOLVE(@(U) HBMRESFUN_TMP([U; W0], MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, tol), U0, nopts);
