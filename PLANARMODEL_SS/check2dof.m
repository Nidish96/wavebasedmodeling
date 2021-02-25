clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/HARMONIC/')

%%
M = [1 0;
     0 1];
K = [1 0;
     0 16];
C = [0.005 0;
     0 0.005];
 
Fv = [1; 0.0];

mu = 0.6;
kt = 3;
kn = 9;
gap = 1;
 
MDL = MDOFGEN(M, K, C, eye(2));

% fnl = @(t, u, varargin) ELDRYFRICT2D(t, u, kt, kn, mu, gap, varargin{:});
% MDL = MDL.SETNLFUN(2+5, eye(2), fnl, eye(2));

%% HBM - Force Controlled
h = [0 1];
Nhc = sum((h==0)+2*(h~=0));

Nt = 2^8;

Fl = zeros(Nhc*2, 1);
Fl(3:4) = Fv;

Wst = 0.1;
Wen = 2.0;

dw = 0.1;
dsmin = 0.0001;
dsmax = 0.5;

Copt = struct('Nmax', 1000, 'itDisplay', false, 'angopt', 1e0, ...
    'dsmin', dsmin, 'arclengthparm', 'arclength', 'dsmax', dsmax, ...
    'lsrch', 0, 'DynDscale', 1);

% Copt.Dscale = [ones(Nhc*2, 1); 1];

U0 = zeros(Nhc*2,1);
UC = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), U0, Wst, Wen, dw, Copt);

% %%
% Sopt = struct('jac', 'full', 'dynamicDscale', 1);
% UC = solve_and_continue(U0, @(Uw) MDL.HBRESFUN(Uw, Fl, h, Nt, 1e-6), Wst, Wen, dw, Sopt);

%
figure(20)
% clf()
npl=semilogy(UC(end,:), sum(UC([3 5],:).^2), 'o-');