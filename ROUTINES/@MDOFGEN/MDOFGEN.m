classdef MDOFGEN
    %MDOFGEN creates a generic mdof model class with given Mass, Stiffness
    %and Damping matrices and specified nonlinearities
    properties
        Ndofs  % Number of DoFs
        
        M   % Mass
        K   % Stiffness
        C   % Damping
        
        L  % displacement transformation matrix
        
        NLTs % Nonlinear functions
        SHAKER % Shaker Attachment
    end
    
    methods
        function m = MDOFGEN(M, K, C, L)
        %MDOFGEN initializes the MDOFGEN object
        % 
        % USAGE:
        %   m = MDOFGEN(M, K, C, L);
            m.M = M;
            m.K = K;
            m.C = C;
            m.L = L;
            
            m.Ndofs = size(m.M,1);
        end
        
        function m = SETNLFUN(m, type, Ls, func, Lf)
        %SETNLFUN sets nonlinear functions to be applied to the model. Two
        %types are allowed: 'inst' and 'hyst' (see below)
        %
        % USAGE:
        %   type    : a+b (possible values: {4, 6, 5, 7}),
        %               with "a" being,
        %     [1] for instantaneous nonlinearity (in time domain) or
        %     [2] for hysteretic nonlinearity (in time domain)
        %               AND "b" being,
        %     [3] for self adjoint force application (F = Ls'*func(Ls*U))
        %     [5] for non-self adjoint force application (F = Lf*func(Ls*u)
        %   Ls      : (Nldofs, Ndofs) selection matrix
        %   func    : if 'inst', [ft, dfdut, dfdudt] = @(t, u, ud) func;
        %                           returning (Nldofs-qtties); assumed to
        %                           be vectorized in time and u (arranged
        %                           properly)
        %             if 'hyst', [ft, dfdut, dfdudt] = @(t, u, tp, up, fp, h)
        %                   func; returning (Nldofs-qtties); assumed to be
        %                   vectorized in u only.
        %   Lf      : [reqd only for type = 7,8] (Ndofs,Nldofs)
        %               "integration" matrix
            nlfun.L    = Ls;
            nlfun.func = func;
            nlfun.type = type;
            if nlfun.type > 5  % Non-self adjoint forcing
                nlfun.Lf = Lf;
            end
            
            m.NLTs = [m.NLTs; nlfun];
        end
        
        function m = ATTACHSHAKER(m, E, A, B, K, Kd, Nshape)
        %ATTACHSHAKER stores the information to attach (one or more) shaker
        %model(s) to the model. Only applicable to explicit (RK) time
        %stepping solvers.
        %  Shaker-attached Model assumed to be written in state-space form as,
        %       E Xd = A X + B U + K*(X-Nshape Y) + Kd*(X - Nshape Yd)
        %       M Ydd + C Yd + K Y + Fnl + Nshape'*K*(Nshape' Y - X) + Nshape'*Kd*(Nshape' Yd - X)
        %   
        %   USAGE:
        %       E, A    : (nX, nX) A matrix 
        %       B       : (nX, nU) B matrix 
        %       K, Kd   : (nX, nX) K matrix 
        %       NshapeK : (nX, nY) NshapeK matrix
            
            m.SHAKER = struct('E', E, 'A', A, 'B', B, 'K', K, 'Kd', Kd, ...
                'Nshape', Nshape);
        end
    end
end

