%% Choices

% exc_pt = 1;
% exc_dir_str = 'y'; % 'x', 'y'
% fsamp = 2^17;
% Wfrc = 259;  % Forcing Frequency (Hz)
% Famp = 1.0;  % Forcing Amplitude

switch exc_dir_str
    case 'x'
        exc_dir = 1;
    case 'y'
        exc_dir = 2;
    otherwise
        error('Unknown direction')
end

%% clc
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
sint = 1e-6;
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

%% Plot Sensor Locations
if false
    %%
    U = zeros(size(Ln,1), 1);
    % U(1:3:end) = 0.01;
    % U(2:3:end) = (1:size(Ln,1)/3)*1e-3;

    figure(1)
    clf()
    set(gcf, 'Color', 'white')
    PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, 1., ccm1, ccm2, alph);
    plot(Xins(:)+Ninvs(1:2:end, :)*U, Ninvs(2:2:end,:)*U, 'ko', 'MarkerFaceColor', 'w')
    xlim([-0.01-0.1 0.73+0.05])
    ylim([-1 1]*0.125)

    % Annotations
    text(Xins(:)-0.01+0.05, zeros(size(Xins(:)))-0.05, sttls)
    for ix=1:length(Xins)
    %     annotation('arrow', (Xins(i)-0.01)*[1 1], [-0.05 0], 'k')
        quiver(Xins(ix)+0.05, -0.04, -0.05, 0.04, 'k', 'LineWidth', 2);
    end

    quiver(BM1.X(1), BM1.Wi, 0, 0.05, 'k', 'ShowArrowHead', 'off')
    quiver(BM2.X(end), BM1.Wi, 0, 0.05, 'k', 'ShowArrowHead', 'off')
    quiver(BM1.X(1)+0.01, BM1.Wi+0.025, BM2.X(end)-BM1.X(1)+0.05, 0, 'k', 'ShowArrowHead', 'off')

    text((BM2.X(end)-BM1.X(1))/2-0.055, BM1.Wi+0.03, sprintf('%d mm', (BM2.X(end)-BM1.X(1))*1000))

    quiver(BM1.X(end), BM1.Wi/2*1.5, 0, BM1.Wi, 'k', 'ShowArrowHead', 'off')
    quiver(BM2.X(1), BM1.Wi/2*1.5, 0, BM1.Wi, 'k', 'ShowArrowHead', 'off')
    quiver(BM1.X(end)+0.01, BM1.Wi*1.1, BM2.X(1)-BM1.X(end)-0.01, 0, 'k', 'ShowArrowHead', 'off')

    text(BM1.X(end)+(BM2.X(1)-BM1.X(end))/2-0.055, BM1.Wi*1.1+0.005, sprintf('%d mm', (BM2.X(1)-BM1.X(end))*1000))

    quiver(BM1.X(1)-0.01, BM1.Wi/2, -0.08, 0, 'k', 'ShowArrowHead', 'off')
    quiver(BM1.X(1)-0.01, -BM1.Wi/2, -0.08, 0, 'k', 'ShowArrowHead', 'off')
    quiver(BM1.X(1)-0.045, -BM1.Wi, 0, BM1.Wi/2, 'k', 'ShowArrowHead', 'off')
    quiver(BM1.X(1)-0.045, BM1.Wi, 0, -BM1.Wi/2, 'k', 'ShowArrowHead', 'off')
    text(BM1.X(1)-0.1, -0.0225, sprintf('%.1f mm', BM1.Wi*1000), 'Rotation', 90)
    axis off
    
    export_fig('./FIGS/BMWSENS_SCHEM.eps', '-depsc')
end

%% Input Force Vector
Finp = Ninvs((exc_pt-1)*2+exc_dir, :)';

%% GMDOF
MDL = MDOFGEN(Ln'*(Mbrb+Mbolt)*Ln, Ln'*(Kbrb+Kbolt)*Ln, Ln'*Cbrb*Ln, Ln);

fnl = @(t, u, varargin) ELDRYFRICT2D(t, u, kt, kn, mu, gap, varargin{:});
MDL = MDL.SETNLFUN(2+5, QrelL, fnl, LTrel);

%% Static Prestress Analysis
U0 = (MDL.K + Ln'*(K0)*Ln)\(Ln'*Fbolt*Prestress);

opts = struct('reletol', 1e-6, 'Display', true);
[Ustat, ~, ~, ~, J0] = NSOLVE(@(U) MDL.STATRESFUN(U, Ln'*Fbolt*Prestress), U0, opts);
U = Ln*Ustat;

if false
    %%
    %U = Lf*Vf(:,3);
    %U = Fbolt;
    sc = 1e3; ccm1 = 'b';ccm2 = 'r';alph = 0.6;

    figure(2)
    clf()
    set(gcf, 'Color', 'white')
    PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph);
%     colorbar('southoutside')
    xlim([-0.01 0.73])
    ylim([-1 1]*0.125)
    xlabel('X Coordinate (m)')
    ylabel('Y Coordinate (m)')
    
    title(sprintf('Displacement Scale: %d x', sc))
    export_fig('./FIGS/PRESTRESSOLN.eps', '-depsc')
end

%% Transient Sine Excitation
% fsamp = 12800;
% fsamp = 100e3;

T0 = 0; T1 = 1.6;
dt = 1/fsamp;
T = (T0:dt:T1)';

%% Excitation
U = zeros(MDL.Ndofs, length(T));
Ud = zeros(MDL.Ndofs, length(T));
Udd = zeros(MDL.Ndofs, length(T));

fprintf('Harmonic Excitation\n F %f N\n Wfrc %f Hz\n Point P%d\n DOF %s\n Sampling %f Hz\n', ...
    Famp, Wfrc, exc_pt, exc_dir, fsamp);
FEX = @(t) Ln'*(Fbolt*Prestress + Finp*Famp*cos(2*pi*Wfrc*t));

%% Initial Condition
E = HARMONICSTIFFNESS(MDL.M, MDL.C, J0, 2*pi*Wfrc, 1);
Ui = E\[Ln'*Finp*Famp; Ln'*Finp*0];

U0  = Ustat + Ui(1:MDL.Ndofs);
Ud0 = (2*pi*Wfrc)*Ui(MDL.Ndofs+(1:MDL.Ndofs));

% U0 = Ustat;
% Ud0 = U0*0;

%% HHTA Nonlinear Implicit Time integration
opts = struct('Display', 'waitbar');
[~, ~, ~, MDL] = MDL.NLFORCE(0, Ustat, zeros(size(Ustat)), 0, 1);

tic 
[T, U, Ud, Udd, MDLh] = MDL.HHTAMARCH(T0, T1, dt, U0, Ud0, ...
                    FEX, opts);
toc

fprintf('===================Done================================\n');

UPs = (Ninvs*Ln)*U;
UdPs = (Ninvs*Ln)*Ud;
UddPs = (Ninvs*Ln)*Udd;
fext = Famp*cos(2*pi*Wfrc*T);

fname = sprintf('./DATS/SSHARM/RESU%d_PT%d%s_F%d_W%d.mat', log2(fsamp), exc_pt, exc_dir_str, Famp*1000, Wfrc);
save(fname, 'T', 'U', 'Ud', 'Udd', 'Finp', 'Famp', 'Wfrc', 'exc_pt', 'exc_dir_str', 'exc_dir', ...
    'Fbolt', 'Prestress', 'UPs', 'UdPs', 'UddPs');

%%
if false
   %%
   figure(4)
   hold on
%    clf()
   plot(T, UddPs((exc_pt-1)*2+exc_dir, :), '.-')
   xlabel('Time (s)')
   ylabel('Displacement at Forcing DOF')
end
