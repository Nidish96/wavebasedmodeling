clc
clear all 
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM/')

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

%% Analytical Wave speeds
% % circular section
% kappa = 10/9;
% nu = 0.29;
% pars = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu))*kappa, ...
%               'A', pi*wdt^2, 'Iz', pi*wdt^4/4);  % Monolithic region 

C1 = pars.rho^2*pars.A*pars.Iz;
C2 = pars.rho*pars.A*pars.Iz*(pars.E+pars.G);
C3 = pars.rho*pars.G*pars.A^2;
C4 = pars.G*pars.E*pars.A*pars.Iz;

% ks = linspace(0, 2*2*pi/(Lmono*2+Lint), 100);
ks = linspace(0, 2*2*pi/wdt, 100);
ws = zeros(4, length(ks));
dwdks = zeros(4, length(ks));
wbyks = zeros(4, length(ks));
for ik=1:length(ks)
  cpol = [C1 0 -C2*ks(ik)^2-C3 0 C4*ks(ik)^4];
  ws(:, ik) = sort(roots(cpol));
  dwdks(:, ik) = (-2*C2*ws(:, ik).^2*ks(ik) + 4*C4*ks(ik)^4)./...
		 (4*C1*ws(:, ik).^3 - 2*C2*ws(:,ik)*ks(ik) - 2*C3*ws(:,ik));
  wbyks(:, ik) = ws(:, ik)/ks(ik);

  if ik==1
    dwdks(~isfinite(dwdks(:,ik)),ik) = 0;
    wbyks(~isfinite(wbyks(:,ik)),ik) = 0;
  end
    
  % fprintf('Done %d\n', ik);
end
c0 = sqrt(pars.E/pars.rho);  % bar velocity

figure(1)
clf();
plot(wdt*ks/2/pi, wbyks(3,:)/c0, 'b-'); hold on 
plot(wdt*ks/2/pi, ks*sqrt(pars.E*pars.Iz/pars.rho/pars.A)/c0, 'k--'); grid on

legend('Timoshenko Beam', 'Euler-Bernoulli Beam')

xlim([0 2]); 
ylim([0 1])
xlabel('ka/2\pi')
ylabel('c/c0')

print('./FIGS/TM_WS.eps', '-depsc')