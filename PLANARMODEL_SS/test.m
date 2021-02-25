clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/PLANARBM/')

%%
wdt = 25.4e-3;
L = 0.72;
Ne = 56;

nu = 0.3;
kappa = 10*(1+nu)/(12+11*nu);
pars = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu))*kappa, ...
              'A', wdt^2, 'Iz', wdt^4/12);  % Monolithic region 

%%          
M = zeros((Ne+1)*3);
K = zeros((Ne+1)*3);
for e=1:Ne
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    [Me, Ke] = PLANTMBMEMATS(L/Ne, pars);
    
    M(is:ie, is:ie) = M(is:ie, is:ie) + Me;
    K(is:ie, is:ie) = K(is:ie, is:ie) + Ke;
end

%%
[V, W] = eig(K, M);
[W, si] = sort(diag(sqrt(W)));
V = V(:, si);
V = V./sqrt(diag(V'*M*V));

disp(W(1:10)/2/pi)

%% Half of a beam
Mb = zeros((Ne/2+1)*3);
Kb = zeros((Ne/2+1)*3);
for e=1:Ne/2
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    [Me, Ke] = PLANTMBMEMATS(L/Ne, pars);
    
    Mb(is:ie, is:ie) = Mb(is:ie, is:ie) + Me;
    Kb(is:ie, is:ie) = Kb(is:ie, is:ie) + Ke;
end

%% Assemble
La = eye((Ne+2)*3);

La((Ne/2)*3+(1:3), (Ne/2+1)*3+(1:3)) = [1 0 -wdt/8;0 1 0;0 0 1];
La(:, (Ne/2)*3+(1:3)) = [];

Ma = La'*blkdiag(Mb, Mb)*La;
Ka = La'*blkdiag(Kb, Kb)*La;

%% 
[Va, Wa] = eig(Ka, Ma);
[Wa, si] = sort(diag(sqrt(Wa)));
Va = Va(:, si);
Va = Va./sqrt(diag(Va'*Ma*Va));

disp(Wa(1:10)/2/pi)
%%
BM = struct('X', linspace(0, L, Ne+1), 'Y', zeros(1, Ne+1), ...
    'Wi', wdt, 'Wo', wdt);

sc = 1e-1;
figure(10)
clf()

plot(BM.X, BM.Y, 'ko-'); hold on
PLANARBMDEPICT(sc*V(:, 4), BM, 'gxy', 0.6)
colorbar('southoutside')

axis equal