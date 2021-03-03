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

exc_pt = 1;
exc_dir_str = 'y';

switch exc_dir_str
    case 'x'
        exc_dir = 1;
    case 'y'
        exc_dir = 2;
    otherwise
        error('Unknown direction')
end

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

% %% Input Force Vector 
% exc_dir_str = 'y';
% Finp = zeros((Nemono+Nein+2)*2*3,1);
% if exc_dir_str=='x'
%     Finp(1) = 1;
% else
%     Finp(2) = 1;
% end
% 
% Finp = [L1;L2]'*Finp;

%% Input/Sensor Locations
X1ins = linspace(BM1.X(1), BM1.X(end), 5); X1ins(end) = [];
N1invs = zeros(length(X1ins)*2, size(Ln,1)); % [ax1; tv1; ax2; tv2; ...]
for ix=1:length(X1ins)
    ein = find((BM1.X(1:end-1)-X1ins(ix)).*(BM1.X(2:end)-X1ins(ix))<=0); 
    ein = ein(1);
    
    V = BM1.X(ein:ein+1);
    xi = -1 + (X1ins(ix)-V(1))*2/(V(2)-V(1));
    
    N1invs((ix-1)*2+1, :) = ([1-xi 1+xi]/2)*L1((ein-1)*3+[1 4], :);
    N1invs((ix-1)*2+2, :) = HERMSF(xi, range(V))*L1((ein-1)*3+[2 3 5 6], :);
end

X2ins = linspace(BM2.X(1), BM2.X(end), 5); X2ins(1) = [];
N2invs = zeros(length(X2ins)*2, size(Ln, 1));
for ix=1:length(X2ins)
    ein = find((BM2.X(1:end-1)-X2ins(ix)).*(BM2.X(2:end)-X2ins(ix))<=0); 
    ein = ein(1);
    
    V = BM2.X(ein:ein+1);
    xi = -1 + (X2ins(ix)-V(1))*2/(V(2)-V(1));
    
    N2invs((ix-1)*2+1, :) = ([1-xi 1+xi]/2)*L2((Nein+1)+(ein-1)*3+[1 4], :);
    N2invs((ix-1)*2+2, :) = HERMSF(xi, range(V))*L2((Nein+1)+(ein-1)*3+[2 3 5 6], :);
end


Xins = [X1ins X2ins];
Ninvs = [N1invs; N2invs];
sttls = cell(length(Xins), 1);
for ix=1:length(Xins)
    sttls{ix} = sprintf('P%d', ix);
end

%% Input Force Vector
Finp = Ninvs((exc_pt-1)*2+exc_dir, :)';

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

if false
    figure(1)
    clf()
    PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph);
    colorbar('southoutside')
    xlim([-0.01 0.73])
    ylim([-1 1]*0.125)
end

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

%% Continuation
Wst = 250*2*pi;
Wen = 270*2*pi;
dw = 1.0;   dsmax = 5;    dsmin = 0.1;

Nt = 2^8;

contin = 0;

Copt = struct('Nmax', 400, 'itDisplay', false, 'angopt', 1e-2, ...
    'dsmin', dsmin, 'adapt', 1, ...
    'arclengthparm', 'orthogonal', 'dsmax', dsmax, ...
    'lsrch', 0, 'DynDscale', 0, 'ITMAX', 100);
% %%
fi = 2;

if contin==0 % Regular stepping
    Copt.adapt = 0;
    Copt.arclengthparm = 'none';
    % Copt.arclengthparm = 'orthogonal';
    dw = 2*pi;
end
UCus = cell(size(Fas));
UCds = cell(size(Fas));
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

    %% UPSWEEP
    [UCus{fi}, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
    if flag==-1
        % Try tangent predictor
        w = UCus{fi}(end,end)+dw;
        [R, dRdU, dRdw] = MDL.HBRESFUN(UCus{fi}(:, end), Fl, h, Nt, 1e-6);
        U0 = UCus{fi}(1:end-1,end) - (dRdU\dRdw)*dw;
        [UC2, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, w, Wen, dw, Copt);
        if flag==-1  % Try linear predictor
            U0 = HARMONICSTIFFNESS(MDL.M, MDL.C, J0, w, h)\Fl;
            U0(1:MDL.Ndofs) = Ustat;
            [UC2, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, w, Wen, dw, Copt);
            if flag==-1 && fi~=1
                previdx = find((UCus{fi-1}(end, 1:end-1)-w).*(UCus{fi-1}(end, 2:end)-w)<=0);
                U0 = UCus{fi-1}(1:end-1, previdx(1));
                U0((MDL.Ndofs+1):end) = U0((MDL.Ndofs+1):end)*(Fas(fi)/Fas(fi-1));
                [UC2, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, w, Wen, dw, Copt);
            end
        end
        
        UCus{fi} = [UCus{fi} UC2];
    end
    % UC = CONTINUE(@(Uw) HBMRESFUN_TMP(Uw, MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

    uout_h = (kron(eye(Nhc), Finp'*Ln)*UCus{fi}(1:end-1,:));
    FRFs{fi}{1} = [(uout_h(2,:)-uout_h(3,:)*1j)/fa; UCus{fi}(end,:)];

    %% DOWNSWEEP
    [UCds{fi}, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wen, Wst, dw, Copt);
    if flag==-1
        % Try tangent predictor
        w = UCds{fi}(end,end)-dw;
        [R, dRdU, dRdw] = MDL.HBRESFUN(UCds{fi}(:, end), Fl, h, Nt, 1e-6);
        U0 = UCds{fi}(1:end-1,end) - (dRdU\dRdw)*(-dw);
        [UC2, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), ...
            U0, w, Wst, dw, Copt);
        if flag==-1  % Try linear predictor
            U0 = HARMONICSTIFFNESS(MDL.M, MDL.C, J0, w, h)\Fl;  % Linear Initial Guess
            U0(1:MDL.Ndofs) = Ustat;
            
            [UC2, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), ...
                U0, w, Wst, dw, Copt);
            if flag==-1 && ~isempty(find((UCus{fi}(end,1:end-1)-w).*(UCus{fi}(end,2:end)-w)<0, 1))  % Try point from upsweep
                U0 = UCus{fi}(1:end-1, (UCus{fi}(end,1:end-1)-w).*(UCus{fi}(end,2:end)-w)<0);
                
                [UC2, ~, ~, flag] = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), ...
                    U0, w, Wst, dw, Copt);
            end
        end
        UCds{fi} = [UCds{fi} UC2];
    end
    % UC = CONTINUE(@(Uw) HBMRESFUN_TMP(Uw, MDL.M, MDL.C, MDL.K, MDL.NLTs.L, MDL.NLTs.Lf, Fl, kt, kn, mu, gap, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

    uout_h = (kron(eye(Nhc), Finp'*Ln)*UCds{fi}(1:end-1,:));
    FRFs{fi}{2} = [(uout_h(2,:)-uout_h(3,:)*1j)/fa; UCds{fi}(end,:)];
    
%    figure(5); plot(FRFs{fi}{1}(2,:)/2/pi, abs(FRFs{fi}{1}(1,:)), '.-', FRFs{fi}{2}(2,:)/2/pi, abs(FRFs{fi}{2}(1,:)), 'o')
end
%%
if contin==0
    save(sprintf('DATS/SSHBM_P%d%s_nocont.mat',exc_pt,exc_dir_str), 'FRFs', 'Fas', 'UCus', 'UCds');
else
    save(sprintf('DATS/SSHBM_P%d%s_withcont.mat',exc_pt,exc_dir_str), 'FRFs', 'Fas', 'UCus', 'UCds');
end
%%
if false
    figure(10)
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
        xlim([Wst Wen]/2/pi)
    
        subplot(2,1, 2)
        % plot(UC(end,:)/2/pi, rad2deg(angle(frf)), '.-'); hold on 
        plot(FRFs{fi}{1}(2,:)/2/pi, rad2deg(angle(FRFs{fi}{1}(1,:))), '.-', 'Color', COLS(fi,:)); hold on
        plot(FRFs{fi}{2}(2,:)/2/pi, rad2deg(angle(FRFs{fi}{2}(1,:))), '.-', 'Color', COLS(fi,:)); hold on
        grid on
        xlabel('Frequency (Hz)')
        ylabel('FRF Phase (degs)')
        xlim([Wst Wen]/2/pi)
    end

    % export_fig('./FIGS/FRF_EXP.eps', '-depsc')
end
