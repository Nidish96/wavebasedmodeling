clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/HARMONIC/')

% set(0,'defaultTextInterpreter','tex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',16)

%% 
fsamp = 2^18;
T0 = 0;
T1 = 1;
T = (0:1/fsamp:T1)';
Nt = length(T);

fmin = 150; 
fmax = 350; 

fp = 250;

ft = sin(2*pi*fp*T);
Np2 = fix(fsamp/fp)*4+2;

wndw = [hanning(Np2); zeros(Nt-Np2, 1)];
% wndw = [ones(Np2,1); zeros(Nt-Np2, 1)];
% wndw = [hamming(Np2); zeros(Nt-Np2, 1)];

[fs, Ffs] = FFTFUN(T, [ft ft.*wndw wndw]);

%% Plotting

figure(1)
clf()
set(gcf, 'Color', 'white')
plot(T, ft, '-', 'LineWidth', 2); hold on;
plot(T, ft.*wndw, '-', 'LineWidth', 2);
plot(T, wndw, 'k-', 'LineWidth', 2)
xlim([0 0.1])

xlabel('Time (s)')
ylabel('Forcing (N)')
export_fig('./FIGS/HNWSIG_T.eps', '-depsc')

figure(2);
clf()
set(gcf, 'Color', 'white')
plot(fp*[1 1], [-300 0], 'LineWidth', 2); hold on
plot(fs, 20*log10(abs(Ffs(:,2))), 'LineWidth', 2)
plot(fs, 20*log10(abs(Ffs(:,3))), 'k-', 'LineWidth', 2)

xlim([0 750])
ylim([-150 0])

xlabel('Frequency (Hz)')
ylabel('Signal (dB)')

legend('Exact sinusoid', 'Windowed Signal', 'Hanning Window')
export_fig('./FIGS/HNWSIG_F.eps', '-depsc')

disp('Done')