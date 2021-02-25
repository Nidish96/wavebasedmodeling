% clc
% clear all 
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/SOLVERS')
addpath('../ROUTINES/HARMONIC')

set(0,'defaultTextInterpreter','latex');
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

[V0, W0] = eig(J0, MDL.M);
[W0, si] = sort(sqrt(diag(W0)));
V0 = V0(:, si);
V0 = V0./sqrt(diag(V0'*MDL.M*V0)');

%% Plot Static Solution
U = Ln*Ustat;

sc = 5e2; ccm1 = 'b';ccm2 = 'r';alph = 0.6;

figure(1)
clf()
set(gcf, 'Color', 'white')

PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph);
axis off

xlim([-0.01 0.73])
% ylim([-0.05 0.245])
ylim([-0.245 0.05])

% ax = axes('Position', [0.225 0.4 0.6 0.55]);
ax = axes('Position', [0.225 0.1 0.6 0.55]);
sc = 5e2; ccm = 'I2';
PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, [1 1 1]*0.8, ccm, 1.0);

colormap(jet)
xx=colorbar('southoutside');
xlabel(xx, 'Strain Second Invariant')
axis off

xlim([IN1.X(1)-0.025 IN1.X(end)+0.025])

export_fig(sprintf('./FIGS/%dIN_PRESSOL.eps', Nein), '-depsc')
disp('Done')

%% Plot Modes
mi = 1; sc = 1e-2; ccm1 = 'b';ccm2 = 'r';alph = 0.6;

figure(2)
clf()
set(gcf, 'Color', 'white')

for mi=1:5
    U = Ln*V0(:, mi);
    
    subplot(5, 1, mi)
    
    if mi==5
        sc = 4e-2
    end
    PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph);
    
    title(sprintf('Mode %d: %.2f Hz', mi, W0(mi)/2/pi), 'fontsize', 20)
    axis equal
    axis off
end

% export_fig(sprintf('./FIGS/%dIN_PRESLINMDS.png', Nein), '-dpng')
disp('Done')

%% Plot on W vs K plot
ks = linspace(0, 2*2*pi*32/120e-3, 1e4); ks(1) = [];
[ws, wbyks, dwdks] = TIMOWSPDS(ks, pars);  % Solution of dispersion relationship

ws = [ws; sqrt(pars.E/pars.rho)*ks];
wbyks = [wbyks; sqrt(pars.E/pars.rho)*ones(size(ks))];

figure(3)
clf()
set(gcf, 'Position', [1170 550 750 420])
set(gcf, 'Color', 'white')

loglog(ks, ws(3:5, :), '-', 'LineWidth', 2); grid on
hold on

% loglog(xlim, sqrt(pars.E/pars.rho)*xlim, 'k--', 'LineWidth', 1.25);
% loglog(xlim, sqrt(pars.G/pars.rho)*xlim, 'k-.', 'LineWidth', 1.25);

symbs = {'-'; '--'; ':'; '-.'};
for iw=1:4
    loglog(xlim, W0(iw)*[1 1], ['k' symbs{iw}])
end
loglog(xlim, W0(5)*[1 1], 'k-', 'LineWidth', 2)

ll = legend('Transverse 1', 'Transverse 2', 'Longitudinal', ...
    'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', 'Location', 'eastoutside');
set(ll, 'Interpreter', 'latex')

xlabel('Wave Number $\kappa=2\pi/\lambda$ (rad/m)')
ylabel('Frequency $\omega=2\pi f$ (rad/s)')

export_fig('./FIGS/TMWS_WvK_wmf.eps', '-depsc')
disp('Done')