clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM/')

% set(0,'defaultTextInterpreter','tex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',16)

Nein = 8;
%% Parameters
% Geometry
Lint = 120e-3;
Lmono = 300e-3;
wdt = 25.4e-3;

nu = 0.3;
kappa = 10*(1+nu)/(12+11*nu);
pars = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu))*kappa, ...
              'A', wdt^2, 'Iz', wdt^4/12);  % Monolithic region 
parsint = pars;  % Interface region 
parsint.A = pars.A/2;
Parsint.Iz = pars.Iz/8;

% Bolt
BoltLocs = [30 60 90]*1e-3;

kappabolt = 6*(1+nu)/(7+6*nu);
rbolt = 8.42e-3/2;
parsbolt = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu))*kappabolt, ...
                  'A', pi*rbolt^2, 'Iz', pi*rbolt^4/4);
[Mbe, Kbe] = PLANTMBMEMATS(25.4e-3, parsbolt);
mhead = 2.864e-2 - parsbolt.rho*parsbolt.A*25.4e-3;
Mbe([1 2 4 5], [1 2 4 5]) = Mbe([1 2 4 5], [1 2 4 5]) + eye(4)*mhead/2;

km = 2.864e-2;  % Member stiffness
Kbe([1 4], [1 4]) = Kbe([1 4], [1 4]) + km*[1 -1;-1 1];

sigbolt = (14.3e-3 + 2*tand(33)*25.4e-3/2)/3;
tbolt = @(x) normpdf(x, BoltLocs(1), sigbolt) + ...
        normpdf(x, BoltLocs(2),sigbolt) + ...
        normpdf(x, BoltLocs(3), sigbolt);  % Traction field

%% FE Discretization
IN1 = struct('X', linspace(Lmono, Lmono+Lint, Nein+1)', ...
             'Y', -wdt/4*ones(Nein+1, 1), ...
             'Wi', wdt/2, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 

IN2 = struct('X', linspace(Lmono, Lmono+Lint, Nein+1)', ...
             'Y', wdt/4*ones(Nein+1, 1), ...
             'Wi', wdt/2, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 

INm = struct('X', linspace(Lmono, Lmono+Lint, Nein+1)', ...
             'Y', zeros(Nein+1, 1), ...
             'Wi', wdt, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 

%% Matrices
M1 = zeros((Nein+1)*3);
K1 = zeros((Nein+1)*3);
M2 = zeros((Nein+1)*3);
K2 = zeros((Nein+1)*3);
Mim = zeros((Nein+1)*3);
Kim = zeros((Nein+1)*3);

for e = 1:Nein
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    
    [Me, Ke] = PLANTMBMEMATS(IN1.X(n2)-IN1.X(n1), parsint);
    
    % IN1
    M1(is:ie, is:ie) = M1(is:ie, is:ie) + Me;
    K1(is:ie, is:ie) = K1(is:ie, is:ie) + Ke;
    
    % IN2
    M2(is:ie, is:ie) = M2(is:ie, is:ie) + Me;
    K2(is:ie, is:ie) = K2(is:ie, is:ie) + Ke;
    
    % IN3
    [Me, Ke] = PLANTMBMEMATS(INm.X(n2)-INm.X(n1), pars);
    Mim(is:ie, is:ie) = Mim(is:ie, is:ie) + Me;
    Kim(is:ie, is:ie) = Kim(is:ie, is:ie) + Ke;
end

%% Bolt Forcing
% Quadrature
Nq = 10;
[~, Tup] = PLANQPMATS(Nq, wdt/4, IN1, parsint);
[~, Tdown, Xqps] = PLANQPMATS(Nq, -wdt/4, IN1, parsint);

tbqps = tbolt(Xqps-IN1.X(1));
scal = 3/sum(Tdown(2:3:end, 2:3:end)'*tbqps);

F1 = scal*(Tdown(2:3:end, :)'*tbqps);
F2 = -scal*(Tup(2:3:end, :)'*tbqps);

figure(1)
clf()
set(gcf, 'Color', 'white')
plot(IN1.X, IN1.X*0, 'ko-'); hold on
plot(Xqps, tbqps, 'b-', 'LineWidth', 2)
for bi=1:3
    text(IN1.X(1)+BoltLocs(bi)-0.008, 41, sprintf('\\textbf{Bolt %d}', bi), 'FontSize', 16)
    plot(IN1.X(1)+BoltLocs(bi)*[1 1], ylim, 'k--')
end

ylabel('Bolt traction ($N/m^2$)')

yyaxis right
stem(IN1.X, abs(F1(2:3:end)), 'k^-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w')

ylabel('Nodal Force (N)')

xlabel('X Coordinate (m)')

ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'k';
% export_fig(sprintf('./FIGS/BOLTFORCMOD_%dIN.eps', Nein), '-depsc')

%% Assembly
Mintf = blkdiag(M1, M2);
Kintf = blkdiag(K1, K2);
Fbolt_intf = [F1; F2];
EYE = eye((Nein+1)*3);
ZER = zeros((Nein+1)*3);

for bi=1:length(BoltLocs)
    % Find element bolt belongs to
    ei = find((IN1.X(1:end-1)-BoltLocs(bi)-IN1.X(1)).*(IN1.X(2:end)-BoltLocs(bi)-IN1.X(1))<=0);
    ei = ei(1);
    
    n1 = ei;
    n2 = ei+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;    
    
    Le = IN1.X(n2)-IN1.X(n1);
    
    x = BoltLocs(bi)+IN1.X(1);
    Nsf = [IN1.X(n2)-x, x-IN1.X(n1)]/Le;
    if abs(max(abs(Nsf))-1)<eps  % Bolt is "at" a node
        Nsf = kron(Nsf, eye(3));
    else % Bolt is in between an element
        Nsf = [Nsf(1), 0, 0, Nsf(2), 0, 0;
               0, Nsf(1), 0, 0, Nsf(2), 0;
               0, 0, 1/2, 0, 0, 1/2];
    end
    
    Qb1 = [1 0 wdt/4;0 1 0;0 0 1]*Nsf*[EYE(is:ie, :) ZER(is:ie, :)];
    Qb2 = [1 0 -wdt/4;0 1 0;0 0 1]*Nsf*[ZER(is:ie, :) EYE(is:ie, :)];
    
    Qbs = [diag([1 -1 1])*Qb1([2 1 3], :); 
           diag([1 -1 1])*Qb2([2 1 3], :)];
    
    Kintf = Kintf + Qbs'*Kbe*Qbs;
    Mintf = Mintf + Qbs'*Mbe*Qbs;
end

%% Quadrature Mats
Nq = 3;
[Qup, Tup] = PLANQPMATS(Nq, wdt/4, IN1, parsint);
[Qdown, Tdown, Xqps] = PLANQPMATS(Nq, -wdt/4, IN1, parsint);

Qrel = [Qup, -Qdown];
Trel = [Tup, -Tdown];

%% Save
save(sprintf('./DATS/INTMATS_%dIN.mat', Nein), 'IN1', 'IN2', 'INm', ...
    'Mintf', 'Kintf', 'Mim', 'Kim', 'Fbolt_intf', 'Qrel', 'Trel', 'Xqps', ...
    'pars', 'parsbolt', 'parsint', 'wdt');

%% Checks
Nemono = fix(Nein*2.5);
BM1 = struct('X', linspace(0, Lmono, Nemono+1)', ...
             'Y', zeros(Nemono+1, 1), ...
             'Wi', wdt, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 
BM2 = struct('X', linspace(Lmono+Lint, Lmono+Lint+Lmono, Nemono+1)', ...
             'Y', zeros(Nemono+1, 1), ...
             'Wi', wdt, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 

Mbm1 = zeros((Nemono+1)*3);
Kbm1 = zeros((Nemono+1)*3);

Mbm2 = zeros((Nemono+1)*3);
Kbm2 = zeros((Nemono+1)*3);

for e=1:Nemono
    % BM1
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    [Me, Ke] = PLANTMBMEMATS(BM1.X(n2)-BM1.X(n1), pars);
    Mbm1(is:ie, is:ie) = Mbm1(is:ie, is:ie) + Me;
    Kbm1(is:ie, is:ie) = Kbm1(is:ie, is:ie) + Ke;
    
    [Me, Ke] = PLANTMBMEMATS(BM2.X(n2)-BM2.X(n1), pars);
    Mbm2(is:ie, is:ie) = Mbm2(is:ie, is:ie) + Me;
    Kbm2(is:ie, is:ie) = Kbm2(is:ie, is:ie) + Ke;    
end

%% Check Monolithic Beam
[Vbm, Wbm] = eig(Kbm1, Mbm1);
[Wbm, si] = sort(sqrt(abs(diag(Wbm))));
Vbm = Vbm(:, si);

kLbm = sqrt(Wbm/sqrt(pars.E*pars.Iz/pars.rho/pars.A))*range(BM1.X);

%%
Mmod = blkdiag(Mbm1, Mbm2);
Kmod = blkdiag(Kbm1, Kbm2);

Lass = eye((Nemono+1)*3*2);
Lass((Nemono+1-1)*3+(1:3), (Nemono+1+1-1)*3+(1:3)) = eye(3);
Lass(:, (Nemono+1-1)*3+(1:3)) = [];

Mmod = Lass'*Mmod*Lass;
Kmod = Lass'*Kmod*Lass;

[V2, W2] = eig(Kmod, Mmod);
[W2, si] = sort(sqrt(abs(diag(W2))));
V2 = V2(:, si);

kL2 = sqrt(W2/sqrt(pars.E*pars.Iz/pars.rho/pars.A))*range(BM1.X)*2;

% Analytical solution (free-free EB Beam): [4.7300 7.8532 10.99564 14.13715 17.2787]

%% Plot em
X2 = BM2.X;
BM2.X = BM2.X - BM2.X(1) + BM1.X(end);

sc = 1e-1; ccm1 = 'b';ccm2 = 'r';alph = 0.6;

mi = 5; 

figure(4); 
clf();
plot(BM1.X, BM1.X*0, 'ko-'); hold on
plot(BM2.X, BM2.X*0, 'ko-'); hold on

PLANARBMDEPICT(Lass(1:(Nemono+1)*3, :)*V2(:, mi)*sc, BM1, ccm1, alph)
PLANARBMDEPICT(Lass((Nemono+1)*3+(1:(Nemono+1)*3), :)*V2(:, mi)*sc, BM2, ccm2, alph)

title(sprintf('Mode %d: %f Hz', mi, W2(mi)/2/pi))

axis equal

BM2.X = X2;
% return

%% Assembly of BRB model
Mmod = blkdiag(Mbm1, Mintf, Mbm2);
Kmod = blkdiag(Kbm1, Kintf, Kbm2);
Fbmod = [zeros((Nemono+1)*3, 1); Fbolt_intf; zeros((Nemono+1)*3, 1)];

Lass = eye((Nemono+1+Nein+1)*3*2);

Mono1_ne = Nemono+1;  % last node of first monolithic region
Mono2_se = Nemono+1+Nein+1+Nein+1+1;  % first node of second monolithic region

Lass((Mono1_ne-1)*3+(1:3), (Mono1_ne+1-1)*3+(1:3)) = [1 0 -wdt/4;
                                                      0 1 0;
                                                      0 0 1];
Lass((Mono2_se-1)*3+(1:3), (Mono2_se-1-1)*3+(1:3)) = [1 0 wdt/4;
                                                      0 1 0;
                                                      0 0 1];

Lass(:, [(Mono1_ne-1)*3+(1:3), (Mono2_se-1)*3+(1:3)]) = [];

Mmod = Lass'*Mmod*Lass;
Kmod = Lass'*Kmod*Lass;
Fbmod = Lass'*Fbmod;

%% 
[V, D] = eig(Kmod, Mmod);
[D, si] = sort(diag(D));
V = V(:, si);

L1 = Lass(1:(Nemono+1+Nein+1)*3, :);
L2 = Lass((Nemono+1+Nein+1)*3+(1:(Nemono+1+Nein+1)*3), :);

sc = 1e-1; ccm1 = 'b';ccm2 = 'r';alph = 0.6;
mi = 4;
U = V(:, mi);

figure(2)
clf()
PLOTSOLN(U, BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph);
colorbar('southoutside')

title(sprintf('Mode %d: %f Hz', mi, sqrt(D(mi))/2/pi))
xlim([-0.01 0.73])
% ylim([-1 1]*0.125)