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
mu   = 0.25;  % 0.25

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

% for i=1:(Nqp*Nein)
%     MDL = MDL.SETNLFUN(2+5, QrelL((i-1)*2+(1:2),:), fnl, LTrel(:, (i-1)*2+(1:2)));
% end

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
h = [0 1];
Nhc = sum(h==0)+2*sum(h~=0);

fa = 25.00;  % 0.01, 0.20, 1.00, 5.00, 15.00, 25.00, 

Fas = [1.0 9.0 17.0 25.0];
FRFs = cell(size(Fas));
%%
Wst = 250*2*pi;
Wen = 270*2*pi;
dw = 1.0;   dsmax = 5;    dsmin = 0.1;

Nt = 2^8;

Copt = struct('Nmax', 400, 'itDisplay', false, 'angopt', 1e-2, ...
    'dsmin', dsmin, 'arclengthparm', 'orthogonal', 'dsmax', dsmax, ...
    'lsrch', 0, 'DynDscale', 0);
% %%
fi = 2;
for fi=1:length(Fas)
    fa = Fas(fi);
    
    Fl = zeros(MDL.Ndofs*Nhc, 1);
    Fl(1:MDL.Ndofs) = Ln'*Fbolt*Prestress;
    Fl(MDL.Ndofs+(1:MDL.Ndofs)) = fa*Ln'*Finp;

    E = HARMONICSTIFFNESS(MDL.M, MDL.C, J0, Wst, h);
    U0 = E\Fl;
    U0(1:MDL.Ndofs) = Ustat;

    % Copt.Dscale = [abs(U0); Wst];
    
    if fi<length(Fas)
        Copt.angopt = 1e-4;
    else
        Copt.angopt = 1e-3;
    end
    
    FRFs{fi} = cell(2,1);

    UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
    % UC = CONTINUE(@(Uw) HBMRESFUN_TMP(Uw, MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

    uout_h = (kron(eye(Nhc), Finp'*Ln)*UC(1:end-1,:));
    FRFs{fi}{1} = [(uout_h(2,:)-uout_h(3,:)*1j)/fa; UC(end,:)];
    
    UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wen, Wst, dw, Copt);
    % UC = CONTINUE(@(Uw) HBMRESFUN_TMP(Uw, MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

    uout_h = (kron(eye(Nhc), Finp'*Ln)*UC(1:end-1,:));
    FRFs{fi}{2} = [(uout_h(2,:)-uout_h(3,:)*1j)/fa; UC(end,:)];    
end
%%
figure(9)
clf(); 
set(gcf, 'Color', 'white')

COLS = DISTINGUISHABLE_COLORS(length(Fas));
for fi=1:length(Fas)
    subplot(2,1, 1)
    % plot(UC(end,:)/2/pi, 20*log10(abs(frf)), '.-'); hold on 
    plot(FRFs{fi}{1}(2,:)/2/pi, 20*log10(abs(FRFs{fi}{1}(1,:))), '.-', 'Color', COLS(fi,:)); hold on
    plot(FRFs{fi}{2}(2,:)/2/pi, 20*log10(abs(FRFs{fi}{2}(1,:))), '.-', 'Color', COLS(fi,:)); hold on
    grid on
    xlabel('Frequency (Hz)')
    ylabel('|FRF| (dB)')

    subplot(2,1, 2)
    % plot(UC(end,:)/2/pi, rad2deg(angle(frf)), '.-'); hold on 
    plot(FRFs{fi}{1}(2,:)/2/pi, rad2deg(angle(FRFs{fi}{1}(1,:))), '.-', 'Color', COLS(fi,:)); hold on
    plot(FRFs{fi}{2}(2,:)/2/pi, rad2deg(angle(FRFs{fi}{2}(1,:))), '.-', 'Color', COLS(fi,:)); hold on
    grid on
    xlabel('Frequency (Hz)')
    ylabel('FRF Phase (degs)')
end

export_fig('./FIGS/FRF_EXP.eps', '-depsc')