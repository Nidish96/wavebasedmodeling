clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
addpath('../ROUTINES/FEM/')

% set(0,'defaultTextInterpreter','tex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',16)

%% Parameters
wdt = 25.4e-3;
Len = 160;
rho = 7800;
E = 2e11;
Ar = wdt^2;

c = sqrt(E/rho)*1.0;

xf = Len/2;

%% FE
Ne = 100;
Xn = linspace(0, Len, Ne+1);
Le = Xn(2);

M = zeros(Ne+1);
K = zeros(Ne+1);
C = zeros(Ne+1);

for e=1:Ne
    M(e:e+1, e:e+1) = M(e:e+1, e:e+1) + rho*Ar*Le/6*[2 1; 1 2];
    K(e:e+1, e:e+1) = K(e:e+1, e:e+1) + Ar*E/Le*[1 -1;-1 1];
    
%     M(e:e+1, e:e+1) = M(e:e+1, e:e+1) + Le/6*[2 1; 1 2];
%     K(e:e+1, e:e+1) = K(e:e+1, e:e+1) + c^2/Le*[1 -1;-1 1];
end
C([1 end],[1 end]) = c*rho*Ar*eye(2);

% Forcing
F = zeros(Ne+1,1);
fe = find((Xn(1:end-1)-xf).*(Xn(2:end)-xf)<=0);
for e=fe
    F(e:e+1) = [Xn(e+1)-xf; xf-Xn(e)]/Le;
end
F = F/sum(F);

%% Null transform
% L = null(ones(Ne+1,1)'*M);
L = eye(size(M));

%% MDOFGEN
MDL = MDOFGEN(L'*M*L, L'*K*L, L'*C*L, L);

%% Transient Analysis
fsamp = 2^12;
T0 = 0;
T1 = 0.5;
T = (T0:1/fsamp:T1)';

famp = 1;
fpuls = 350;
fpwdt = (2*pi)/(2*pi*fpuls/sqrt(E/rho));
ft = sin(2*pi*fpuls*T)*famp;
Np2 = ceil(fsamp/fpuls)*4+1;
% wndw = [zeros(Np2,1); hanning(Np2); zeros(length(T)-Np2*2,1)];
wndw = [hanning(Np2); zeros(length(T)-Np2,1)];

opts = struct('Display', 'waitbar');
U0 = zeros(MDL.Ndofs,1);
Ud0 = zeros(MDL.Ndofs,1);

tic
[T, U, Ud, Udd] = MDL.HHTAMARCH(T0, T1, 1/fsamp, ...
    U0, Ud0, ((L'*F).*(ft(:).*wndw(:))'), ...
    opts);
toc


%% Plot
figure(10)
clf()

plot(T, (L'*F)'*U, '.-')

xlabel('Time (s)')
ylabel('Response')
pause(1)

%%
figure(11)
clf()
pause(0.01)
um = max(max(abs(U)));
for ti=1:length(T)
    clf()
    plot(Xn, Xn*0, 'k.-'); hold on
    
    plot(Xn, L*U(:, ti), '.-');
    
    xlabel('X Coordinate (m)')
    ylabel('Response')
    ylim(um*[-1 1])
    
    title(sprintf('Frame %d/%d. F = %.2f N', ti, length(T), ft(ti)*wndw(ti)))
    
    pause(0.001)
end