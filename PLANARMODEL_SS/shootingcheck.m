



opts = struct('Display', 'none');
% [T, U0, Ud0, Udd0, ~, dUU, dUUd, dUdU, dUdUd, dUddU, dUddUd] = MDL.HHTAMARCH(T0, dt, dt, Ustat, zeros(size(Ustat)), FEX, opts);
dt = 2e-5;

T0 = 0;
T1 = 1/Wfrc;
dt = T1/fix(T1/dt);

UUd0 = [Ustat; zeros(size(Ustat))];
tic
[R, dRdUUd] = SHRESFUN(UUd0, MDL, T0, T1, dt, FEX, opts);
toc

% it = 0;
% while true
%     [R, dRdUUd] = SHRESFUN(UUd0, MDL, T0, T1, dt, FEX, opts);
%     UUd1 = UUd0 - dRdUUd\R;
%     
%     it = it+1;
%     fprintf('%d %e %e %e\n', it, vecnorm(R), vecnorm(UUd1-UUd0), (UUd1-UUd0)'*R)
%     
%     UUd0 = UUd1;
% end


fopts = optimset('Jacobian', 'off', 'Display', 'iter');  % options for fsolve if NSOLVE fails 
UUdsol = fsolve(@(uud) SHRESFUN(uud, MDL, T0, T1, dt, FEX, opts), UUd0, fopts);

%%
function [R, dRdUUd] = SHRESFUN(UUd, MDL, T0, T1, dt, FEX, opts)
    [T, U, Ud, Udd, ~, dUU, dUUd, dUdU, dUdUd, dUddU, dUddUd] = ...
        MDL.HHTAMARCH(T0, T1, dt, UUd(1:MDL.Ndofs), UUd(MDL.Ndofs+(1:MDL.Ndofs)), FEX, opts);
    
    R = [U(:,end)-UUd(1:MDL.Ndofs);
        Ud(:,end)-UUd(MDL.Ndofs+(1:MDL.Ndofs))];
    dRdUUd = [dUU-eye(MDL.Ndofs), dUUd;
        dUdU, dUdUd-eye(MDL.Ndofs)];
end