clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC')
addpath('../ROUTINES/SOLVERS')

%% Description
% Using a 2dof model to determine if the implemented functions for ELDRYFRICT2D provide correct gradients
% dof 1 -> tangential; dof 2 -> normal

%% Model
M = [1 0;0 1];
K = [4 0;0 4];
C = [2*0.001*2 0;0 2*0.001*2];
    
mu = 0.9;
kt = 5;
kn = 5;

MDL = MDOFGEN(M, K, C, eye(2));

fnl = @(t, u, varargin) ELDRYFRICT2D(t, u, kt, kn, mu, 0, varargin{:});
MDL = MDL.SETNLFUN(2+5, eye(2), fnl, eye(2));

%% Static Analysis
Fn = [0; 1];
Ft = [1; 0];

% Static forces
fns = 10;
fts = 0;

U0 = MDL.K\(Fn*fns+Ft*fts);

opts = struct('reletol', 1e-6, 'Display', true);
[Ustat, ~, ~, ~, J0] = NSOLVE(@(U) MDL.STATRESFUN(U, Fn*fns+Ft*fts), U0, opts);

%% HBM
h = [0 1 2 3];
Nhc = sum((h==0)+2*(h~=0));

Nd = 2;
Nt = 2^7;

Wst = 0.1;
Wen = 5;

dw = 0.1;

% Dynamic forces
fnd = 0;
ftd = 1;

Fl = [Fn*fns+Ft*fts;
    Fn*fnd+Ft*ftd;
    zeros(Nd*(Nhc-2),1)];

U0 = HARMONICSTIFFNESS(MDL.M,MDL.C,J0,Wst,h)\Fl;

Copt = struct('Nmax', 100, 'itDisplay', false, 'angopt', 5e-1,...
    'arclengthparm', 'arclength', 'lsrch', 2);
% Copt.Dscale = [ones(Nd*Nhc,1)*5e0; 1];
UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

% Nw = 400;
% UC = zeros(Nd*Nhc+1, Nw);
% UC(end,:) = linspace(Wst, Wen, Nw);
% opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'off');
% u0 = HARMONICSTIFFNESS(MDL.M,MDL.C,J0,UC(end,1),h)\Fl;
% 
% for iw=1:Nw
%     [UC(1:end-1,iw), ~, eflg, ~, Jac] = fsolve(@(u) MDL.HBRESFUN([u; UC(end,iw)], Fl, h, Nt, 1e-6), u0, opts);
%     
%     if eflg<=0
%         keyboard
%     end
%     
%     [~,Jac,Jw] = MDL.HBRESFUN(UC(:,iw), Fl, h, Nt, 1e-6);
%     if iw<Nw
%         u0 = UC(1:end-1,iw) - (Jac\Jw)*(UC(end,iw+1)-UC(end,iw));
%     end
%     fprintf('%d/%d: %d\n', iw,Nw, eflg);
% end

figure(100)
clf();
subplot(2,1, 1)
plot(UC(end,:), abs(UC(3:4,:)+UC(5:6,:)*1j), '.-')
xlabel('Frequency (Hz)')
ylabel('Response (m)')

subplot(2,1, 2)
plot(UC(end,:), rad2deg(angle(UC(3:4,:)+UC(5:6,:)*1j)), '.-')
xlabel('Frequency (Hz)')
ylabel('Phase (degs)')

%% 
iw = 57;

opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'iter', 'CheckGradients', false, 'MaxIterations', 1000);
[~,Jac,Jw] = MDL.HBRESFUN(UC(:,iw-1), Fl, h, Nt, 1e-6);
u0 = UC(1:end-1, iw-1) - (Jac\Jw)*(UC(end,iw)-UC(end,iw-1));
% u0 = HARMONICSTIFFNESS(MDL.M,MDL.C,J0,UC(end,iw),h)\Fl;
% u0 = UC(1:end-1,iw-1);
UC(1:end-1,iw) = fsolve(@(u) MDL.HBRESFUN([u; UC(end,iw)], Fl, h, Nt, 1e-6), u0, opts);

%% Check gradient
[R0, dR0, dRw0] = MDL.HBRESFUN([u0; UC(end,iw)], Fl, h, Nt, 1e-6);
dRnum = zeros(size(dR0));
dRwnum = zeros(size(R0));

hv = zeros(Nd*Nhc, 1);
hm = 1e-6;

for hi = 1:Nd*Nhc
    hv(hi) = 1;
    Rp = MDL.HBRESFUN([u0+hv*hm; UC(end,iw)], Fl, h, Nt, 1e-6);
    Rm = MDL.HBRESFUN([u0-hv*hm; UC(end,iw)], Fl, h, Nt, 1e-6);
    hv(hi) = 0;
    
    dRnum(:, hi) = (Rp-Rm)/(2*hm);
end
Rp = MDL.HBRESFUN([u0; UC(end,iw)+hm], Fl, h, Nt, 1e-6);
Rm = MDL.HBRESFUN([u0; UC(end,iw)-hm], Fl, h, Nt, 1e-6);
dRwnum = (Rp-Rm)/(2*hm);

disp('dR0')
disp(dR0)
disp('dRnum')
disp(dRnum)

disp('dRw0')
disp(dRw0)
disp('dRwnum')
disp(dRwnum)