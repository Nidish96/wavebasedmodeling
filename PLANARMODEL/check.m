clc
clear all 
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxesFontSize',16)

%% Parameters
% Geometry
Lint = 120e-3;
Lmono = 300e-3;
wdt = 25.4e-3;

BoltLocs = [30 60 90]*1e-3;

nu = 0.3;
kappa = 10*(1+nu)/(12+11*nu);
pars = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu))*kappa, ...
              'A', wdt^2, 'Iz', wdt^4/12);  % Monolithic region 
parsint = pars;  % Interface region 
parsint.A = pars.A/2;
Parsint.Iz = pars.Iz/8;

%% Finite Element Discretization
Nein = 8;  % 8 elements in interface 
Nemono = fix(Nein*2.5);   % elements in monolithic parts 

BM1 = struct('X', linspace(0, Lmono, Nemono+1)', ...
             'Y', zeros(Nemono+1, 1), ...
             'Wi', wdt, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 
IN1 = struct('X', linspace(Lmono, Lmono+Lint, Nein+1)', ...
             'Y', -wdt/4*ones(Nein+1, 1), ...
             'Wi', wdt/2, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 

IN2 = struct('X', linspace(Lmono, Lmono+Lint, Nein+1)', ...
             'Y', wdt/4*ones(Nein+1, 1), ...
             'Wi', wdt/2, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 
BM2 = struct('X', linspace(Lmono+Lint, Lmono+Lint+Lmono, Nemono+1)', ...
             'Y', zeros(Nemono+1, 1), ...
             'Wi', wdt, ... % in-plane width 
             'Wo', wdt);  % out-of-plane width 

##figure(1)
##clf()
##plot(BM1.X, BM1.Y, 'ko-'); hold on 
##plot(IN1.X, IN1.Y, 'bo-'); 
##plot(BM2.X, BM2.Y, 'kx-'); hold on 
##plot(IN2.X, IN2.Y, 'rx-')
##
##grid on 

Mm = zeros((Nemono+1)*3);
Km = zeros((Nemono+1)*3);

for e=1:Nemono
  n1 = e;
  n2 = e+1;
  
  is = (n1-1)*3+1;
  ie = n2*3;
  
  [Me, Ke] = PLANTMBMEMATS(BM1.X(n2)-BM1.X(n1), pars);
  Mm(is:ie, is:ie) = Mm(is:ie, is:ie) + Me;
  Km(is:ie, is:ie) = Km(is:ie, is:ie) + Ke;
end

Lb = eye((Nemono+1)*3);
% Lb(:, [1 2 3]) = [];  % fixed free 
Lb(:, [1 2 end-2 end-1]) = [];   % pinned pinned 

Mm = Lb'*Mm*Lb;
Km = Lb'*Km*Lb;

[Vs, Ws] = eig(Km, Mm);
[Ws, si] = sort(sqrt(diag(Ws)));

Vs = Lb*Vs(:, si);

mi = 1;
sc = 4e-2;
figure(1)
clf();
plot(BM1.X, BM1.X*0, 'ko-'); hold on; 
grid on 

% PLANARBMDEPICT(Vs(:, mi)*sc, BM1, 'b', 0.1);

PLANARBMDEPICT(Vs(:, mi)*sc, BM1, 'ex', 0.8);
colorbar('southoutside')

xlim([0 0.35])
ylim([-1 1]*0.35/2)