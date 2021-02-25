clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')

c = 5;
%% 
xl0 = [sqrt(c); 0];

l0 = 0;
l1 = 6;

ds = 0.01;

Copt = struct('Nmax', 200, 'angopt', 1e2, ...
    'arclengthparm', 'arclength');
UC = CONTINUE(@(xl) resfun(xl, c), xl0(1), l0, l1, ds, Copt);

figure(1)
clf()
plot(UC(end,:), UC(1,:), '.-'); hold on
xas = linspace(sqrt(c), -sqrt(c), 100);
las = c-xas.^2;
plot(las, xas, 'k--')

%%
function [r, drdx, drdl] = resfun(xl, c)    
    r = xl(1)^2+xl(2)-c;
    drdx = 2*xl(1);
    drdl = 1;
end