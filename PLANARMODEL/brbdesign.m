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

%##figure(1)
%##clf()
%##plot(BM1.X, BM1.Y, 'ko-'); hold on 
%##plot(IN1.X, IN1.Y, 'bo-'); 
%##plot(BM2.X, BM2.Y, 'kx-'); hold on 
%##plot(IN2.X, IN2.Y, 'rx-')
%##
%##grid on 

%axis equal 
xlim([-0.1 0.4])
ylim([-0.25 0.25])


%% Matrices 
M1 = zeros((Nemono+Nein+2)*3);
M2 = zeros((Nemono+Nein+2)*3);
K1 = zeros((Nemono+Nein+2)*3);
K2 = zeros((Nemono+Nein+2)*3);

for e=1:Nemono
    % BM1
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    [Me, Ke] = PLANTMBMEMATS(BM1.X(n2)-BM1.X(n1), pars);
    M1(is:ie, is:ie) = M1(is:ie, is:ie) + Me;
    K1(is:ie, is:ie) = K1(is:ie, is:ie) + Ke;
    
    % BM2 
    n1 = Nein+1+e;
    n2 = Nein+1+e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    [Me, Ke] = PLANTMBMEMATS(BM2.X(e+1)-BM2.X(e), pars);
    M2(is:ie, is:ie) = M2(is:ie, is:ie) + Me;
    K2(is:ie, is:ie) = K2(is:ie, is:ie) + Ke;
end

for e=1:Nein
  % IN1 
  n1 = Nemono+1+e;
  n2 = Nemono+1+e+1;
  
  is = (n1-1)*3+1;
  ie = n2*3;
  [Me, Ke] = PLANTMBMEMATS(IN1.X(e+1)-IN1.X(e), parsint);
  M1(is:ie, is:ie) = M1(is:ie, is:ie) + Me;
  K1(is:ie, is:ie) = K1(is:ie, is:ie) + Ke;
  
  % IN2
  n1 = e;
  n2 = e+1;
  
  is = (n1-1)*3+1;
  ie = n2*3;
  [Me, Ke] = PLANTMBMEMATS(IN2.X(n2)-IN2.X(n1), parsint);
  M2(is:ie, is:ie) = M2(is:ie, is:ie) + Me;
  K2(is:ie, is:ie) = K2(is:ie, is:ie) + Ke;
end 

%% Assembly
L1 = eye((Nemono+Nein+2)*3);
indms1 = Nemono*3+(1:3);
indis1 = (Nemono+1)*3+(1:3);
L1(indms1, indis1) = [1 0 -wdt/4;0 1 0;0 0 1];
L1(:, indms1) = [];

M1t = L1'*M1*L1;
K1t = L1'*K1*L1;

L2 = eye((Nemono+Nein+2)*3);
indis2 = Nein*3+(1:3);
indms2 = (Nein+1)*3+(1:3);
L2(indms2, indis2) = [1 0 wdt/4;0 1 0;0 0 1];
L2(:, indms2) = [];

M2t = L2'*M2*L2;
K2t = L2'*K2*L2; 

%% Assemble together
Mbrb = blkdiag(M1t, M2t);
Kbrb = blkdiag(K1t, K2t);
L1 = sparse([L1 zeros(size(L2))]);
L2 = sparse([zeros(size(L2)) L2]);

%% Nodal displacements (and rotations)
Ri1 = zeros((Nein+1)*3, size(M1,1));
Ri1(:, (Nemono+1)*3+(1:(Nein+1)*3))      = kron(eye(Nein+1), [1 0 -wdt/4;0 1 0;0 0 1]);
Ri2 = zeros((Nein+1)*3, size(M1,1));
Ri2(:, 1:(Nein+1)*3)      = kron(eye(Nein+1), [1 0 wdt/4;0 1 0;0 0 1]);

%% Relative displacements at Quadrature points
Nqp = 10;
Q1 = zeros((Nein*Nqp)*2, (Nein+1)*3);
Q2 = zeros((Nein*Nqp)*2, (Nein+1)*3);
T1 = zeros((Nein*Nqp)*2, (Nein+1)*3);
T2 = zeros((Nein*Nqp)*2, (Nein+1)*3);
Xqps = zeros(Nein*Nqp, 1);
for e=1:Nein
  n1 = e;
  n2 = e+1;
  
  is = (n1-1)*3+1;
  ie = n2*3;
  
  Le = IN1.X(n2)-IN1.X(n1);
  
  [xi, wi] = LGWT(Nqp, 0, Le); 
  xi = xi(end:-1:1);
  wi = wi(end:-1:1);
  Xqps((e-1)*Nqp+(1:Nqp)) = IN1.X(e)+xi;
  
  wi = wi*IN1.Wo;  % Multiplying with out of plane width 
  for qi=1:Nqp
    Q1((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 -wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2];
    Q2((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2];
    
    T1((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 -wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2]*wi(qi);
    T2((e-1)*(Nqp*2)+(qi-1)*2+(1:2), is:ie) = [1 0 wdt/4;0 1 0]*[1-xi(qi)/Le 0 0 xi(qi)/Le 0 0;0 1-xi(qi)/Le 0 0 xi(qi)/Le 0;0 0 1/2 0 0 1/2]*wi(qi);
  end 
end

Qrel = Q1*L1((Nemono+1)*3+(1:(Nein+1)*3), :) - Q2*L2(1:(Nein+1)*3, :);
Trel = (T1*L1((Nemono+1)*3+(1:(Nein+1)*3), :) - T2*L2(1:(Nein+1)*3, :))';

%% Bolt Modeling 
kappa = 6*(1+nu)/(7+6*nu);
rad = 8.42e-3/2;
parsbolt = struct('E', 2e11, 'rho', 7800, 'G', 2e11/(2*(1+nu))*kappa, ...
              'A', pi*rad^2, 'Iz', pi*rad^4/4);
[Mbe, Kbe] = PLANTMBMEMATS(25.4e-3, parsbolt);
mhead = 2.864e-2 - sum(diag(Mbe(1:3:end,1:3:end)))*6/4;
Mbe([1 2 4 5], [1 2 4 5]) = Mbe([1 2 4 5], [1 2 4 5]) + eye(4)*mhead;

km = 2.864e-2;  % Member stiffness in the normal direction
Kbe([1 4], [1 4]) = Kbe([1 4], [1 4]) + km*[1 -1;-1 1];

Mbolt = zeros(size(Mbrb));
Kbolt = zeros(size(Kbrb));
Fbolt = zeros(size(Kbrb,1),1);
for bi=1:length(BoltLocs)
  ei = find((IN1.X(1:end-1)-BoltLocs(bi)-IN1.X(1)).*(IN1.X(2:end)-BoltLocs(bi)-IN1.X(1))<=0);
  ei = ei(1);
  
  n1 = ei;
  n2 = ei+1;
  Le = IN1.X(n2)-IN1.X(n1);
  
  is = (n1-1)*3+1;
  ie = n2*3;
  
  x = BoltLocs(bi)+IN1.X(1);
  Nsf = [IN1.X(n2)-x, x-IN1.X(n1)]/Le;
  
  Qb1 = [1 0 wdt/4;0 1 0;0 0 1]*kron(Nsf, eye(3))*L1((Nemono+1)*3+((ei-1)*3+(1:6)), :);
  Qb2 = [1 0 -wdt/4;0 1 0;0 0 1]*kron(Nsf, eye(3))*L2((ei-1)*3+(1:6), :);
  
  Qbs = [diag([1 -1 1])*Qb1([2 1 3], :); diag([1 -1 1])*Qb2([2 1 3], :)];
  
  Kbolt = Kbolt + Qbs'*Kbe*Qbs;
  Mbolt = Mbolt + Qbs'*Mbe*Qbs;
  
  %% Bolt forcing 
  Fbolt = Fbolt + Qb1(2,:)'-Qb2(2,:)';
end 

%% Bolt forcing through top traction distribution (gaussian)
sig = (14.3e-3 + 2*tand(33)*25.4e-3/2)/3;
tbolt = @(x) normpdf(x, BoltLocs(1), sig)+normpdf(x, BoltLocs(2), sig)+normpdf(x, BoltLocs(3), sig);
tbqps = tbolt(Xqps-IN1.X(1));
scal = 3./sum((T2(2:2:end,:)*L1((Nemono+1)*3+(1:(Nein+1)*3), :))'*tbqps);  % Normalize to unit force from each bolt 
Fboltd = scal*(T2(2:2:end,:)*L1((Nemono+1)*3+(1:(Nein+1)*3), :)-T1(2:2:end,:)*L2(1:(Nein+1)*3, :))'*tbqps;

save('./DATS/BRBMATS.mat', 'Mbrb', 'Kbrb', 'Kbolt', 'Mbolt', 'Fbolt', 'Qrel', ...
            'Trel', 'BM1', 'IN1', 'BM2', 'IN2', 'Nein', 'Nemono', 'pars', ...
            'parsint', 'Ri1', 'Ri2', 'L1', 'L2', 'wdt', 'Nqp', 'Q1', 'T1', ...
            'Q2', 'T2', 'BoltLocs', 'Lmono', 'Lint', 'nu', 'Xqps', 'Fboltd', '-v7');

% Fixed Interface
%##Rm = Ri1*L1 - Ri2*L2; 
%##Lfix = null(Rm([1:3:end 2:3:end],:));  

% Stiff compliant interface 
krel = zeros(Nein*Nqp*2,1);
krel(1:2:end) = 1e6;  % tangential
krel(2:2:end) = 1e9;  % normal 
Kstiff = Trel*diag(krel)*Qrel;
%#Lfix = eye(size(Kbrb));  
%#Kbrb = Kbrb+Kstiff;

% Bolted Interface
Lfix = eye(size(Kbrb));
Kbrb = Kbrb+Kbolt;
Mbrb = Mbrb+Mbolt;

K = Lfix'*Kbrb*Lfix;
M = Lfix'*Mbrb*Lfix; 

[Vf, Ws] = eig(K, M);
[Ws, si] = sort(sqrt(diag(Ws)));
Vf = Vf(:, si);

Vs = Lfix*Vf;

mi = 4;
sc = 1e-1;
figure(1)
clf();
plot(BM1.X, BM1.Y, 'ko-'); hold on 
plot(IN1.X, IN1.Y, 'ko-'); hold on
plot(BM2.X, BM2.Y, 'ko-'); hold on 
plot(IN2.X, IN2.Y, 'ko-'); hold on 
grid on 

%#PLANARBMDEPICT(L1(1:(Nemono+1)*3, :)*Vs(:, mi)*sc, BM1, 'b', 0.1);
%#PLANARBMDEPICT(L1((Nemono+1)*3+(1:(Nein+1)*3), :)*Vs(:, mi)*sc, IN1, 'b', 0.1);
%#
%#PLANARBMDEPICT(L2((Nein+1)*3+(1:(Nemono+1)*3), :)*Vs(:, mi)*sc, BM2, 'r', 0.1);
%#PLANARBMDEPICT(L2(1:(Nein+1)*3, :)*Vs(:, mi)*sc, IN2, 'r', 0.1);

PLANARBMDEPICT(L1(1:(Nemono+1)*3, :)*Vs(:, mi)*sc, BM1, 'gxy', 0.6);
PLANARBMDEPICT(L1((Nemono+1)*3+(1:(Nein+1)*3), :)*Vs(:, mi)*sc, IN1, 'gxy', 0.6);

PLANARBMDEPICT(L2((Nein+1)*3+(1:(Nemono+1)*3), :)*Vs(:, mi)*sc, BM2, 'gxy', 0.6);
PLANARBMDEPICT(L2(1:(Nein+1)*3, :)*Vs(:, mi)*sc, IN2, 'gxy', 0.6);
colorbar('southoutside')

xlim([-0.01 0.73])
ylim([-1 1]*0.125)

%% 
Lnull = null(Vf(:, 1:3)'*M);
Bs = Lnull*((Lnull'*(Kbrb+Kbolt+Kstiff)*Lnull)\(Lnull'*Fbolt*12e3));

%Bs = Fbolt;
sc = 1e2
ccm = 'p1';
figure(2)
clf();
plot(BM1.X, BM1.Y, 'ko-'); hold on 
plot(IN1.X, IN1.Y, 'ko-'); hold on
plot(BM2.X, BM2.Y, 'ko-'); hold on 
plot(IN2.X, IN2.Y, 'ko-'); hold on 
grid on 

PLANARBMDEPICT(L1(1:(Nemono+1)*3, :)*Bs*sc, BM1, ccm, 0.6);
PLANARBMDEPICT(L1((Nemono+1)*3+(1:(Nein+1)*3), :)*Bs*sc, IN1, ccm, 0.6);

PLANARBMDEPICT(L2((Nein+1)*3+(1:(Nemono+1)*3), :)*Bs*sc, BM2, ccm, 0.6);
PLANARBMDEPICT(L2(1:(Nein+1)*3, :)*Bs*sc, IN2, ccm, 0.6);
colorbar('southoutside')

xlim([-0.01 0.73])
ylim([-1 1]*0.125)
