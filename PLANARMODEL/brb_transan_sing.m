%clc
% clear all 
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/SOLVERS')
addpath('../ROUTINES/HARMONIC')

set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxesFontSize',16)

load('./DATS/BRBMATS.mat', 'Mbrb', 'Kbrb', 'Kbolt', 'Mbolt', 'Fbolt', 'Fboltd', 'Qrel', ...
            'Trel', 'BM1', 'IN1', 'BM2', 'IN2', 'Nein', 'Nemono', 'pars', ...
            'parsint', 'Ri1', 'Ri2', 'L1', 'L2', 'wdt', 'Nqp', 'Q1', 'T1', ...
            'Q2', 'T2', 'BoltLocs', 'Lint', 'Lmono', 'nu')
Prestress = 12e3;
Fbolt = Fboltd;

%% Getting rid of fixed interface modes 
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
sint = 5e-6;
chi  = 2;
ktkn = chi*(1-nu)/(2-nu);
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn   = kt/ktkn;

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
Finp = zeros((Nemono+Nein+2)*2*3,3);
Finp(1:3, :) = eye(3);

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

%U = Lf*Vf(:,3);
%U = Fbolt;
sc = 1e3; ccm1 = 'b';ccm2 = 'r';alph = 0.6;

figure(1)
clf()
PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph);
colorbar('southoutside')
xlim([-0.01 0.73])
ylim([-1 1]*0.125)

%% Transient impulse response 
% fsamp = 12800;
T0 = 0; T1 = 0.05;
dt = 1/fsamp;
T = (T0:dt:T1)';

%% Excitation
switch DOF
    case 'X'
        Fin = Finp(:, 1);
    case 'Y'
        Fin = Finp(:, 2);
    case 'T'
        Fin = Finp(:, 3);
end

% Haversine impulse 
type = 'IMP';
% bw = 1000;
% famp = 1000;
% fex = @(t) famp*sin(2*pi*bw*t).^2.*(t<=1.0/(2*bw));
% fext = fex(T);
% FEX = @(t) Ln'*Finp*fex(t)+Ln'*Fbolt*Prestress;

% White Gaussian Noise
type = 'WGN';
% famp = 10;  % 0.01
bw = -1;

rand("seed", sd);

fprintf('Simulating with type %s, F %d, exciting %s, sampled at %d Hz\n', type, famp, DOF, famp);

fext = wgn(length(T), 1, 40+20*log10(famp))';
FEX = @(t) Ln'*Fin*interp1(T, fext, t)+Ln'*Fbolt*Prestress;

%% HHTA Nonlinear Implicit Time integration
opts = struct('Display', 'waitbar');

FEXv = (Ln'*Fin).*fext + Ln'*Fbolt*Prestress;

[~, ~, ~, MDL] = MDL.NLFORCE(0, Ustat, zeros(size(Ustat)), 0, 1);

tic 
[T, U, Ud, Udd, MDL] = MDL.HHTAMARCH(T0, T1, dt, Ustat, zeros(size(Ustat)), ...
                    FEXv, opts);
toc

fname = sprintf('./DATS/%dIN_%sRESP_%s%d_samp%d_r%d.mat', Nein, type, DOF, famp, log2(fsamp), sd);
% save(fname, 'T', 'U', 'Ud', 'Udd', 'fext', 'Finp', 'Fin', 'Nrep');

%%
Wlin = sqrt(sort(eig(J0, MDL.M)))/2/pi;

figure(2)
clf()
plot(T, (Ln'*Finp)'*Udd, '.-')
% plotyy(T, (Ln'*Finp)'*Udd, T, fext)

xlabel('Time (s)')
ylabel('Response')

udd = (Ln'*Finp)'*Udd;
[freqs, Uf] = FFTFUN(T, udd');
[freqs, Ff] = FFTFUN(T, fext');

figure(3)
clf()
subplot(2,1, 1)
plot(freqs, 20*log10(abs(Uf./Ff))); hold on 
grid on
% for i=1:10
%   plot(Wlin(i)*[1 1], ylim, 'k--')
% end 
ylabel('|FRF| (dB)')

subplot(2,1, 2)
plot(freqs, rad2deg(angle(Uf./Ff))); hold on 
% for i=1:10
%   plot(Wlin(i)*[1 1], ylim, 'k--')
% end 
grid on

xlabel('Frequency (Hz)')
ylabel('Phase (degs)')