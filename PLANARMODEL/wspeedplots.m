clc
clear all

set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',16)


load('./DATS/INTMATS_8IN.mat', 'pars');
Le = 120e-3/32;

%% Calculation
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

loglog(ks, sqrt(pars.E/pars.rho).*ks', 'k--', 'LineWidth', 1.25);
loglog(ks, sqrt(pars.G/pars.rho).*ks', 'k-.', 'LineWidth', 1.25);

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