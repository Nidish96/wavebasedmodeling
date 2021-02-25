function [outs, wbyks, dwdks] = TIMOWSPDS(ins, pars, varargin)
%TIMOWSPDS returns the solution of the dispersion relationship for all "k"
%
%  USAGE:
%     [ws, wbyks, dwdks] = TIMOWSPDS(ks, pars);
%  INPUTS:
%       ins     : (N, 1) Spatial/temporal frequencies
%       pars    : structure with E, rho, G, A, Iz
%       which   : 'k'[default] or 'w' (which is input)
%  OUTPUTS:
%       outs    : (4, N) Temporal/spatial frequencies
%       wbyks   : (4, N) w/k
%       dwdks   : (4, N) dw/dk
        
    which = 'k';
    if nargin>=3
        which = varargin{1};
    end

    if which=='k'
        ks = ins;
        
        C1 = pars.rho^2*pars.A*pars.Iz;
        C2 = pars.rho*pars.A*pars.Iz*(pars.E+pars.G);
        C3 = pars.rho*pars.G*pars.A^2;
        C4 = pars.G*pars.E*pars.A*pars.Iz;

        ws = zeros(4, length(ks));
        wbyks = zeros(4, length(ks));
        dwdks = zeros(4, length(ks));
        for ik=1:length(ks)
            cpol = [C1 0 -C2*ks(ik)^2-C3 0 C4*ks(ik)^4];
            ws(:, ik) = sort(roots(cpol));
            dwdks(:, ik) = (-2*C2*ws(:, ik).^2*ks(ik) + 4*C4*ks(ik)^4)./...
             (4*C1*ws(:, ik).^3 - 2*C2*ws(:,ik)*ks(ik) - 2*C3*ws(:,ik));
            wbyks(:, ik) = ws(:, ik)/ks(ik);

            if ik==1
                dwdks(~isfinite(dwdks(:,ik)),ik) = 0;
                wbyks(~isfinite(wbyks(:,ik)),ik) = 0;
            end
        end
        
        outs = ws;
    elseif which=='w'        
        ws = ins; 
        
        C1 = pars.G*pars.E*pars.A*pars.Iz;
        C2 = pars.rho*pars.A*pars.Iz*(pars.E+pars.G);
        C3 = pars.rho^2*pars.A*pars.Iz;
        C4 = pars.rho*pars.G*pars.A^2;

        ks = zeros(4, length(ws));
        wbyks = zeros(4, length(ws));
        dkdws = zeros(4, length(ws));
        for ik=1:length(ws)
            cpol = [C1 0 C2*ws(ik)^2 0 C3*ws(ik)^4-C4*ws(ik)^2];
            ks(:, ik) = sort(roots(cpol));
            dkdws(:, ik) = (C4*ws(ik)-2*C3*ws(ik)^3)./(C2*ks(:, ik)+2*C1*ks(:, ik).^3);
            wbyks(:, ik) = ws(ik)./ks(ik);

            if ik==1
                dkdws(~isfinite(dkdws(:,ik)),ik) = 0;
                wbyks(~isfinite(wbyks(:,ik)),ik) = 0;
            end
        end
        
        outs = ks;
        dwdks = 1./dkdws;
    else
        error('unknown');
    end
end