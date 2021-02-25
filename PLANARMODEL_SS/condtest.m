clc
clear all

%%
nopts = struct('reletol', 1e-6, 'rtol', 1e-6, 'etol', 1e-6, 'utol', 1e-6, ...
    'Display', true, 'lsrch', 0, 'crit', 30)
% nopts.Dscale = [1;1e-7];
x0 = single([0; 0]);
[Xs,~,~,~,~,reu_wo] = NSOLVE(@(x) RESFUN(x), x0, nopts);

nopts.Dscale = [1;1e-7];
[Xs,~,~,~,~,reu_ws] = NSOLVE(@(x) RESFUN(x), [0; 0], nopts);

%%
i = 1;
figure(1)
clf()
semilogy(sqrt(reu_wo(:, i)), 'o-'); hold on
semilogy(sqrt(reu_ws(:, i)), 'o-'); hold on

%% 
function [R, dR] = RESFUN(x)
    sc = single(1e12);
    
    R = single([(x(1)-1)^2;
        (sc*x(2)-0.75)^2]);        
    dR = single([2*(x(1)-1) 0;
        0 2*(sc*x(2)-0.75)*sc]);
end