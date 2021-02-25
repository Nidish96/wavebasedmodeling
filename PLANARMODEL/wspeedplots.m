clc
clear all

set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',16)

%% Calculation
load('./DATS/INTMATS_8IN.mat', 'pars');
Le = 120e-3/32;
ks = linspace(0, 2*2*pi/Le, 1e4); ks(1) = [];
[ws, wbyks, dwdks] = TIMOWSPDS(ks, pars);  % Solution of dispersion relationship

ws = [ws; sqrt(pars.E/pars.rho)*ks];
wbyks = [wbyks; sqrt(pars.E/pars.rho)*ones(size(ks))];


%% Plotting
figure(6)
clf()
set(gcf, 'Color', 'white')

loglog(ks, ws(3:5, :), '-', 'LineWidth', 2); grid on
hold on

% loglog(ks, sqrt(pars.E/pars.rho).*ks', 'k--', 'LineWidth', 1.25);
% loglog(ks, sqrt(pars.G/pars.rho).*ks', 'k-.', 'LineWidth', 1.25);
loglog(xlim, sqrt(pars.E/pars.rho)*xlim, 'k--', 'LineWidth', 1.25);
loglog(xlim, sqrt(pars.G/pars.rho)*xlim, 'k-.', 'LineWidth', 1.25);

% ll = legend('Transverse 1', 'Transverse 2', 'Longitudinal', ...
%     '$\sqrt{E/\rho}$', '$\sqrt{G/\rho}$', 'Location', 'southeast');
% set(ll, 'Interpreter', 'latex')

xlabel('Wave Number $\kappa=2\pi/\lambda$ (rad/m)')
ylabel('Frequency $\omega=2\pi f$ (rad/s)')
export_fig('./FIGS/TMWS_WvK.eps', '-depsc')

figure(7)
clf()
set(gcf, 'Color', 'white')

loglog(ks, wbyks(3:5, :), '-', 'LineWidth', 2); grid on
hold on
loglog(ks, sqrt(pars.E/pars.rho).*ones(size(ks))', 'k--', 'LineWidth', 1.25);
loglog(ks, sqrt(pars.G/pars.rho).*ones(size(ks))', 'k-.', 'LineWidth', 1.25);

ll = legend('Transverse 1', 'Transverse 2', 'Longitudinal', ...
    '$\sqrt{E/\rho}$', '$\sqrt{G/\rho}$', 'Location', 'southeast');
set(ll, 'Interpreter', 'latex')

xlabel('Wave Number $\kappa=2\pi/\lambda$ (rad/m)')
ylabel('Wave Speed $c=\omega/\kappa$ (m/s)')
export_fig('./FIGS/TMWS_CvK.eps', '-depsc')

figure(8)
clf()
set(gcf, 'Color', 'white')

for i=3:5
    loglog(ws(i,:), wbyks(i, :), '-', 'LineWidth', 2); hold on
end
grid on
loglog(xlim, sqrt(pars.E/pars.rho)*[1 1], 'k--', 'LineWidth', 1.25);
loglog(xlim, sqrt(pars.G/pars.rho)*[1 1], 'k-.', 'LineWidth', 1.25);

% ll = legend('Transverse 1', 'Transverse 2', 'Longitudinal', ...
%     '$\sqrt{E/\rho}$', '$\sqrt{G/\rho}$', 'Location', 'southeast');
% set(ll, 'Interpreter', 'latex')

xlabel('Frequency $\omega=2\pi f$ (rad/s)')
ylabel('Wave Speed $c=\omega/\kappa$ (m/s)')
export_fig('./FIGS/TMWS_CvW.eps', '-depsc')

disp('Done')

%% Calculations of necessary domain size

load('./DATS/INTMATS_8IN.mat', 'pars');
Le = 120e-3/32;
ks = linspace(0, 2*2*pi/Le, 1e4);
[ws, wbyks, dwdks] = TIMOWSPDS(ks, pars);  % Solution of dispersion relationship

ws = [ws; sqrt(pars.E/pars.rho)*ks];
wbyks = [wbyks; sqrt(pars.E/pars.rho)*ones(size(ks))];


fmin = 150;
fmax = 350;

k1s = interp1(ws(3, :), ks, 2*pi*[fmin fmax]);  % Transverse 1
k2s = interp1(ws(4, :), ks, 2*pi*[fmin fmax]);  % Transverse 2
k3s = interp1(ws(5, :), ks, 2*pi*[fmin fmax]);  % Longitudinal

c1s = interp1(ws(3, :), wbyks(3, :), 2*pi*[fmin fmax]);  % Transverse 1
c2s = interp1(ws(4, :), wbyks(4, :), 2*pi*[fmin fmax]);  % Transverse 2
c3s = interp1(ws(5, :), wbyks(5, :), 2*pi*[fmin fmax]);  % Longitudinal

gs1s = interp1(ws(3, :), dwdks(3, :), 2*pi*[fmin fmax]);  % Transverse 1
gs2s = interp1(ws(4, :), dwdks(4, :), 2*pi*[fmin fmax]);  % Transverse 2
gs3s = interp1(ws(5, :), wbyks(5, :), 2*pi*[fmin fmax]);  % Longitudinal

fprintf('Frequency (Hz) & Wavenumber ($rad m^{-1}$) & Wavelength (m) & Phase Speed (m/s) & Group Speed (m/s)\\\\\n')
fprintf('%.2f & (%.2f,%.2f,%.2f) &  (%.2f,%.2f,%.2f) &  (%.2f,%.2f,%.2f) &  (%.2f,%.2f,%.2f)\\\\\n', ...
    fmin, k1s(1), k2s(1), k3s(1), 2*pi/k1s(1), 2*pi/k2s(1), 2*pi/k3s(1), c1s(1), c2s(1), c3s(1), gs1s(1), gs2s(1), gs3s(1));
fprintf('%.2f & (%.2f,%.2f,%.2f) &  (%.2f,%.2f,%.2f) &  (%.2f,%.2f,%.2f) &  (%.2f,%.2f,%.2f)\\\\\n', ...
    fmax, k1s(2), k2s(2), k3s(2), 2*pi/k1s(2), 2*pi/k2s(2), 2*pi/k3s(2), c1s(2), c2s(2), c3s(2), gs1s(2), gs2s(2), gs3s(2));