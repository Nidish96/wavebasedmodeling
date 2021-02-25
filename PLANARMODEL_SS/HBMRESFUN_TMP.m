function [R, dRdU, dRdw] = HBMRESFUN_TMP(Uw, M, C, K, Qxn, Txn, Fl, kt, kn, mu, gap, h, Nt, tol)

    Nd = size(M, 1);
    Nhc = sum((h==0)+2*(h~=0));
    
    w = Uw(end);
    
    t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
    cst = zeros(Nt, Nhc);
    if h(1)==0
        cst(:, [1 2:2:end]) = cos(h(:)'.*t);
        cst(:, 3:2:end) = sin(h(2:end)'.*t);
    else
        cst(:, 1:2:end) = cos(h(:)'.*t);
        cst(:, 1:2:end) = sin(h(:)'.*t);
    end
        

%% Linear Part
    [E, dEdw] = HARMONICSTIFFNESS(M, C, K, w, h);
    
%% Nonlinear Friction Forces
    Uxn = Qxn*reshape(Uw(1:end-1), Nd, Nhc);  % [u10 u1ai u1bi; u20 u2ai u2bi; ...]
    uxnt = TIMESERIES_DERIV(Nt, h, Uxn', 0);  % [ut(t) ut(2) ...]
    
    Nnl = size(Qxn, 1)/2;
    
    if length(kt)<=Nnl || length(kn)<=Nnl || length(mu)<=Nnl || length(gap)<=Nnl
        kt  = kt(1)*ones(Nnl,1);
        kn  = kn(1)*ones(Nnl,1);
        mu  = mu(1)*ones(Nnl,1);
        gap = gap(1)*ones(Nnl,1);
    end
    
    fn  = zeros(Nt, Nnl);
    ft  = zeros(Nt, Nnl);
    jnn = zeros(Nt, Nnl, Nnl);
    jtt = zeros(Nt, Nnl, Nnl*Nhc);
    jtn = zeros(Nt, Nnl, Nnl*Nhc);
    
    FNL = zeros(Nnl*Nhc*2, 1);
    JNL = zeros(Nnl*Nhc*2, Nnl*Nhc*2);
    for di=1:Nnl
        xi = (di-1)*2+1;
        ni = (di-1)*2+2;
        
        fn(:, di) = kn(di)*max(uxnt(:, ni)-gap(di), 0);
        jnn(:, di, di) = kn(di)*(fn(:, di)~=0);
        
        its = 0;
        fprev = ft(end, di); 
        
        while (abs(ft(end, di)-fprev)>tol) || (its==0)
            fprev = ft(end, di); 
            
            for ti=1:Nt
                tim1 = mod(ti-1-1, Nt)+1;

                if fn(ti, di)==0  % Separation
                    ft(ti, di) = 0;
                    jtt(ti, di, di:Nnl:end) = 0;
                    jtn(ti, di, di:Nnl:end) = 0;
                else
                    ftpred = kt(di)*(uxnt(ti, xi)-uxnt(tim1, xi))+ft(tim1, di);
                    fslip = mu(di)*fn(ti, di);

                    if abs(ftpred)<=fslip  % Stick
                        ft(ti, di) = ftpred;

                        jtt(ti, di, di:Nnl:end) = kt(di)*(cst(ti, :)-cst(tim1, :))' + squeeze(jtt(tim1, di, di:Nnl:end));
                        jtn(ti, di, di:Nnl:end) = jtn(tim1, di, di:Nnl:end);
                    else  % Slip
                        ft(ti, di) = sign(ftpred)*fslip;

                        jtt(ti, di, di:Nnl:end) = 0;
                        jtn(ti, di, di:Nnl:end) = sign(ftpred)*mu(di)*kn(di)*cst(ti, :);
                    end
                end
            end
            its = its+1;
        end
        
        FNL(xi:Nnl*2:end) = GETFOURIERCOEFF(h, ft(:, di));
        FNL(ni:Nnl*2:end) = GETFOURIERCOEFF(h, fn(:, di));
        
        JNL(ni:Nnl*2:end, ni:Nnl*2:end) = GETFOURIERCOEFF(h, jnn(:, di, di).*cst);
        
        JNL(xi:Nnl*2:end, xi:Nnl*2:end) = GETFOURIERCOEFF(h, squeeze(jtt(:, di, di:Nnl:end)));
        JNL(xi:Nnl*2:end, ni:Nnl*2:end) = GETFOURIERCOEFF(h, squeeze(jtn(:, di, di:Nnl:end)));
    end
    
%% Residual Assembly
    R = E*Uw(1:end-1) + kron(eye(Nhc), Txn)*FNL - Fl;
    dRdU = E + kron(eye(Nhc), Txn)*JNL*kron(eye(Nhc), Qxn);
    dRdw = dEdw*Uw(1:end-1);
    
%     dRdU = [dRdU dRdw];
    
end