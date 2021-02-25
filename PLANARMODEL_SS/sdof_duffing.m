clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/TRANSIENT/')
addpath('../ROUTINES/HARMONIC/')

%% Parameters
m = 1.0;
k = 4.0;
c = 2.0/(2*sqrt(k*m));
b = 0.10;

MDL = MDOFGEN(m, k, c, 1.0);
fnl = @(t,u,ud) deal(b*u.^3, 3*b*u.^2, zeros(size(u)));
MDL = MDL.SETNLFUN(1+3, 1, fnl);

%% Force Controlled HBM
Nt = 256;
h = [0 1 2 3 4 5];  Nhc = sum(h==0)+2*sum(h~=0);
Fl = [0 0 1 zeros(1, Nhc-3)]';
U0 = zeros(Nhc,1);

fa = 40;

Wst = 0.01;
Wen = 5.0;

dw = 0.02;

E = HARMONICSTIFFNESS(m, c, k, Wst, h);
U0 = E\Fl;
Dscale = max([abs(U0); Wst], 0.1);
Copt = struct('Nmax', 1000, 'Dscale', Dscale, 'DynDscale', 1, 'angopt', 1e-1, ...
    'arclengthparm', 'arclength');

% Dscale = ones(size(Dscale));
% Copt.Dscale = Dscale;   
% Copt.DynDscale = 0;
UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl*fa, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

figure(1)
% clf()
hold on
plot(UC(end,:), sqrt([1 0.5*ones(1,Nhc-1)]*UC(1:end-1,:).^2), '.-')

xlabel('Frequency')
ylabel('Amplitude')

% %% Amplitude Controlled HBM
% Nt = 256;
% h = [0 1 3];  Nhc = sum(h==0)+2*sum(h~=0);
% Fl = [0 0 1 zeros(1, Nhc-3)]';
% U0 = zeros(Nhc,1);
% 
% Amp = 10;
% 
% Wst = 0.1;
% Wen = 4.0;
% 
% dw = 0.25;
% 
% E = HARMONICSTIFFNESS(m, c, k, Wst, h);
% U0 = E\Fl;
% Dscale = max([abs(U0); 1; Wst], 0.1);
% Copt = struct('Nmax', 1000, 'Dscale', Dscale, 'DynDscale', 1, 'angopt', 1e-1);
% % UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl*fa, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);
% 
% cons = struct('order', 1, 'type', 'D', 'mask', 1, 'sc', 1e0);
% UC = CONTINUE(@(Ufw) MDL.ACHBRESFUN(Ufw, Amp, Fl, cons, h, Nt, 1e-6), [U0; 0.01], Wst, Wen, dw, Copt);
% 
% figure(2)
% % clf()
% hold on
% plot(UC(end,:), sqrt([1 0.5*ones(1,Nhc-1)]*UC(1:end-2,:).^2)./UC(end-1,:), '.-')
% % plot(UC(end,:), sqrt([1 0.5*ones(1,Nhc-1)]*UC(1:end-2,:).^2), '.-')
% % plot(UC(end,:), sqrt(sum(UC(2:3,:).^2)), '.-')
% % plot(UC(end,:), UC(end-1,:), '.-')
% 
% xlabel('Frequency')
% ylabel('Amplitude')
% 
% figure(3)
% % clf()
% hold on
% plot3(UC(end,:), sqrt(sum(UC(2:3,:).^2)), abs(UC(end-2,:)), '.-')
% grid on
% 
% xlabel('Frequency')
% ylabel('Amplitude |X|')
% zlabel('Force |F|')