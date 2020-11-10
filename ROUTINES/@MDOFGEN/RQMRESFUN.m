function [R, dRdUl, dRdq] = RQMRESFUN(m, Ulq, lsc)

    if lsc==1
        Q = 10^Ulq(end);
        dQdlq = Q*log(10);
    else
        Q = Ulq(end);
        dQdlq = 1.0;
    end
    
    FNL = zeros(m.Ndofs,1);
    dFNL = zeros(m.Ndofs);
    for ni=1:length(m.NLTs)
        if mod(m.NLTs(ni).type,2)==0
            [Fnl, dFnl, ~] = m.NLTs(ni).func(0,  m.NLTs(ni).L*Ulq(1:end-2));  % Additional arguments ignored, implying zeros
        else
            [Fnl, dFnl] = m.NLTs(ni).func(0,  m.NLTs(ni).L*Ulq(1:end-2));  % Additional arguments ignored, implying zeros
        end
        
        if m.NLTs(ni).type<=5  % Self-adjoint forcing
            FNL = FNL + m.NLTs(ni).L'*Fnl;
            dFNL = dFNL + m.NLTs(ni).L'*dFnl*m.NLTs(ni).L;
        else
            FNL = FNL + m.NLTs(ni).Lf*Fnl;
            dFNL = dFNL + m.NLTs(ni).Lf*dFnl*m.NLTs(ni).L;
        end
    end
    
%     % Residue - "Vanilla" version
%     R = [m.K*Ulq(1:end-2)+FNL-Ulq(end-1)*m.M*Ulq(1:end-2);
%         0.5*Ulq(1:end-2)'*m.M*Ulq(1:end-2)-0.5*Q^2];
%     dRdUl = [m.K+dFNL-Ulq(end-1)*m.M, -m.M*Ulq(1:end-2);
%         Ulq(1:end-2)'*m.M, 0];
%     dRdq = [zeros(m.Ndofs,1);-Q*dQdlq]; 
    
    % Residue - Better conditioned version (Jacobian not nearly singular at
    % solution)
    R = [m.K*Ulq(1:end-2)+FNL-Ulq(end-1)*m.M*Ulq(1:end-2);
        Ulq(1:end-2)'*(m.K*Ulq(1:end-2)+FNL)-Ulq(end-1)*Q^2];
    
    dRdUl = [m.K+dFNL-Ulq(end-1)*m.M, -m.M*Ulq(1:end-2);
        Ulq(1:end-2)'*(2*m.K+dFNL)+FNL', -Q^2];
    dRdq = [zeros(m.Ndofs,1);-2*Ulq(end-1)*Q*dQdlq]; 
end