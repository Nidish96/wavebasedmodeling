function [T, U, Ud, Udd, m, varargout] = HHTAMARCH(m, T0, T1, dt, U0, Ud0, Fex, varargin)
%HHTAMARCH conducts time series marching using HHTA. Member of
%MDOFGEN class.
%
%  USAGE:
%    [T, U, Ud, Udd, m] = m.HHTAMARCH(T0, T1, dt, U0, Ud0, Fex, opts);
%       (or)
%    [T, U, Ud, Udd, m, dUU, dUUd, dUdU, dUdUd, dUddU, dUddUd] = ...
%                             m.HHTAMARCH(T0, T1, dt, U0, Ud0, Fex, opts);
%  INPUTS:
%    T0, T1, dt	: Starting, ending, and interval times
%    U0, Ud0 	: Ndx1 displacement and velocity initial conditions
%    Fex	: (1x1)->(Ndx1) forcing function
%    opts	: [optional] Structure with parameters
%    	alpha [0], beta [1/4], gamma [1/2], reletol [1e-6], etol
%    	[1e-6], utol [1e-6], rtol [1e-6], Display ['progress'],
%    	ITMAX [10]
%  OUTPUTS:
%    T		: 1xNt time vector
%    U, Ud, Udd	: NdxNt displacement, velocity, and acceleration 
%    m		: MDOFGEN struct (with updated hysteretic states if
%                 relevant)

% [average acceleration]: alpha[0], beta[1/4], gamma[1/2]
% [implicit linear acceleration]: alpha[0], beta[1/6], gamma[1/2]
% stability condition: 2b >= g >= 1/2
% recom: alpha \in [-1/3, 0]; beta = (1-alpha^2)/4; gamma = 1/2-alpha

  opts = struct('alpha', 0, 'beta', 1/4, 'gamma', 1/2, ... % Newmark
      'reletol', 1e-6, 'etol', 1e-6, 'utol', 1e-6, 'rtol', 1e-6, ...
      'Display', 'progress', ...  % can be 'iter', 'progress', 'both', 'waitbar'
      'ITMAX', 10);
  if length(varargin)==1
      nflds = fieldnames(varargin{1});
      for i = 1:length(nflds)
        opts.(nflds{i}) = varargin{1}.(nflds{i});
      end
  end
  
%   opts.alpha = -1/6;
%   opts.beta  = (1-opts.alpha)^2/4;
%   opts.gamma = 1/2-opts.alpha;
  
  a = opts.alpha;
  b = opts.beta;
  g = opts.gamma;
  
  
  U0  = reshape(U0, m.Ndofs, 1);
  Ud0 = reshape(Ud0, m.Ndofs, 1);
  
  Z1 = m.M + (1+a)*g*dt*m.C + (1+a)*b*dt^2*m.K;
  Z2 = m.M - (1+a)*(1-g)*dt*m.C - (1+a)*(0.5-b)*dt^2*m.K;
  Z3 = (1+a)*dt*m.K;
  
  T = (T0:dt:T1)';    Nt = length(T);
  U = zeros(m.Ndofs, Nt);  Ud = U;  Udd = U;
  U(:, 1) = U0;  Ud(:, 1) = Ud0;
  
  % Initialize acceleration
  [Fnl, dFnldu, dFnldud, m] = m.NLFORCE(T0, U0, Ud0, T0-dt);
%   if strcmp(typeinfo(Fex), 'function handle')
%     Udd(:,1) = m.M\(Fex(T0)-m.C*Ud0-m.K*U0-Fnl);
%   elseif strcmp(typeinfo(Fex), 'matrix')
%     Udd(:,1) = m.M\(Fex(:,1)-m.C*Ud0-m.K*U0-Fnl);
%   end 
  
  if strcmp(class(Fex), 'function_handle')
    Udd(:,1) = m.M\(Fex(T0)-m.C*Ud0-m.K*U0-Fnl);
  elseif strcmp(class(Fex), 'double') || strcmp(class(Fex), 'single')
    Udd(:,1) = m.M\(Fex(:,1)-m.C*Ud0-m.K*U0-Fnl);
  else
      error('wth');
  end
  
  if nargout==11
    dUim1_Udd0 = -(m.K+dFnldu)\m.M;
    dUdim1_Udd0 = -(m.C+dFnldud)\m.M;
    
    dUim1_U0 = eye(m.Ndofs);
    dUim1_Ud0 = zeros(m.Ndofs);
    
    dUdim1_U0 = zeros(m.Ndofs);
    dUdim1_Ud0 = eye(m.Ndofs);
    
    dUddim1_U0 = m.M\(-m.K-dFnldu);
    dUddim1_Ud0 = m.M\(-m.C-dFnldud);
    
    dFnldu0 = dFnldu;
    dFnldud0 = dFnldud;
  end
  clear U0 Ud0
  
  if strcmp(opts.Display, 'waitbar')
      wb = waitbar(T(1)/T1, sprintf('Progress: %.4e/%.4e dt: %.4e', T(1), T1, dt), ...
          'createcancelbtn', "setappdata(gcbf, 'interrupt', true)");
  end
  for i=2:Nt
      %% Explicit Predictor
      Udd(:, i) = Udd(:, i-1);
      
      %% Corrector Iterations
      [FnlP, dFnldu, dFnldud, ~] = m.NLFORCE(T(i-1)+(1+a)*dt, ...
          U(:, i-1) + (1+a)*dt*Ud(:, i-1) + (1+a)*dt^2*((.5-b)*Udd(:, i-1)+b*Udd(:,i)), ...
          Ud(:, i-1) + (1+a)*dt*((1-g)*Udd(:, i-1)+g*Udd(:, i)), T(i-1));
      % Residual, Jacobian, and Update
%       if strcmp(typeinfo(Fex), 'function handle')
%         R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
%             (FnlP-Fnl) - (Fex(T(i-1)+(1+a)*dt)-Fex(T(i-1)));  
%       elseif strcmp(typeinfo(Fex), 'matrix')
%         R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
%             (FnlP-Fnl) - (1+a)*(Fex(:,i)-Fex(:,i-1));
%       else
%         error('wth');        
%       end
      if strcmp(class(Fex), 'function_handle')
        R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
            (FnlP-Fnl) - (Fex(T(i-1)+(1+a)*dt)-Fex(T(i-1)));  
      elseif strcmp(class(Fex), 'double') || strcmp(class(Fex), 'single')
        R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
            (FnlP-Fnl) - (1+a)*(Fex(:,i)-Fex(:,i-1));
      else
        error('wth');
      end

      J = Z1 + (1+a)*(b*dt^2*dFnldu + g*dt*dFnldud);

      if ~isempty(m.NLTs)  % Nonlinear Case, Need to iterate
          du = -J\R;
          % Error norms
          e = abs(R'*du);
          r = mean(R.^2);
          u = mean(du.^2);

          e0 = e;
          r0 = r;
          u0 = u;
          it = 0;          
          
          flag = 8*(e/e0<opts.reletol) + 4*(e<opts.etol) + 2*(r<opts.rtol) + ...
              1*(u<opts.utol);
          if strcmp(opts.Display, 'iter') || strcmp(opts.Display, 'both')
              fprintf('ITN, E, E/E0, r, du\n%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, flag);
              %fprintf('---------------------------------------------------\n');
          end
          while ((flag<7) || (it==0))
              Udd(:, i) = Udd(:, i) + du;
              it = it+1;

              [FnlP, dFnldu, dFnldud, ~] = m.NLFORCE(T(i-1)+(1+a)*dt, ...
                  U(:, i-1) + (1+a)*dt*Ud(:, i-1) + (1+a)*dt^2*((.5-b)*Udd(:, i-1)+b*Udd(:,i)), ...
                  Ud(:, i-1) + (1+a)*dt*((1-g)*Udd(:, i-1)+g*Udd(:, i)), T(i-1));
              % Residual, Jacobian, and Updates
              if strcmp(class(Fex), 'function_handle')
                R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
                    (FnlP-Fnl) - (Fex(T(i-1)+(1+a)*dt)-Fex(T(i-1)));
              elseif strcmp(class(Fex), 'double') || strcmp(class(Fex), 'single')
                R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
                    (FnlP-Fnl) - (1+a)*(Fex(:,i)-Fex(:,i-1));
              else
                error('wth');
              end
              J = Z1 + (1+a)*(b*dt^2*dFnldu + g*dt*dFnldud);
              du = -J\R;
              % Error norms
              e = abs(R'*du);
              r = mean(R.^2);
              u = mean(du.^2);

              flag = 8*(e/e0<opts.reletol) + 4*(e<opts.etol) + 2*(r<opts.rtol) + ...
                  1*(u<opts.utol);
              if strcmp(opts.Display, 'iter') || strcmp(opts.Display, 'both')
                  fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, flag);
              end

              if it>opts.ITMAX
                  flag = 0;
                  break;
              end
          end
          if strcmp(opts.Display, 'iter') || strcmp(opts.Display, 'both')
              fprintf('---------------------------------------------------\n');
          end      
      else  % Linear Case, No need iterations
          flag = 8+4+2+1;
          if strcmp(class(Fex), 'function_handle')
              Udd(:, i) = Z1\(Z2*Udd(:, i-1) - Z3*Ud(:, i-1) - ...
                  (FnlP-Fnl) + (Fex(T(i-1)+(1+a)*dt)-Fex(T(i-1))));
              
%               R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
%                   (FnlP-Fnl) - (Fex(T(i-1)+(1+a)*dt)-Fex(T(i-1)));
          elseif strcmp(class(Fex), 'double') || strcmp(class(Fex), 'single')
              Udd(:, i) = Z1\(Z2*Udd(:, i-1) - Z3*Ud(:, i-1) - ...
                  (FnlP-Fnl) + (1+a)*(Fex(:,i)-Fex(:,i-1)));
              
%               R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
%                   (FnlP-Fnl) - (1+a)*(Fex(:,i)-Fex(:,i-1));
          else
              error('wth');
          end
      end
      
      if flag == 0 || any(~isfinite(abs(Udd(:, i))))
          fprintf('No Convergence/Non finite march at %f s : Returning\n', T(i))
%          keyboard
          
          U = U(:, 1:i-1);
          Ud = Ud(:, 1:i-1);
          Udd = Udd(:, 1:i-1);
          T = T(1:i-1);
          
          if nargout==11
            varargout{1} = dUi_U0;
            varargout{2} = dUi_Ud0;
            varargout{3} = dUdi_U0;
            varargout{4} = dUdi_Ud0;
            varargout{5} = dUddi_U0;
            varargout{6} = dUddi_Ud0;
          end
          
          break;
      end
      
      %% Update States  
      % (alpha used)
      Ud(:, i) = Ud(:, i-1) + (1+a)*dt*((1-g)*Udd(:, i-1)+g*Udd(:, i));
      U(:, i) = U(:, i-1) + (1+a)*dt*Ud(:, i-1) + (1+a)*dt^2*((.5-b)*Udd(:, i-1)+b*Udd(:,i));
      
      % (alpha not used here)
%       Ud(:, i) = Ud(:, i-1) + dt*((1-g)*Udd(:, i-1)+g*Udd(:, i));
%       U(:, i) = U(:, i-1) + dt*Ud(:, i-1) + dt^2*((0.5-b)*Udd(:, i-1)+b*Udd(:, i));
      
      [Fnl, dFnldu, dFnldud, m] = m.NLFORCE(T(i), U(:, i), Ud(:, i), T(i-1));
      
      if nargout==11
%         R = Z1*Udd(:, i) - Z2*Udd(:, i-1) + Z3*Ud(:, i-1) + ...
%             (FnlP-Fnl) - (Fex(T(i-1)+(1+a)*dt)-Fex(T(i-1)));          
        dRdU0 = -Z2*dUddim1_U0 + Z3*dUdim1_U0 + ...
            (dFnldu*(dUim1_U0 + (1+a)*dt*dUdim1_U0 + (1+a)*dt^2*(.5-b)*dUddim1_U0) +...
            dFnldud*(dUddim1_U0 + (1+a)*dt*(1-g)*dUddim1_U0) -...
            dFnldu0*dUim1_U0 - dFnldud0*dUdim1_U0);
        
        dRdUd0 = -Z2*dUddim1_Ud0 + Z3*dUdim1_Ud0 + ...
            (dFnldu*(dUim1_Ud0 + (1+a)*dt*dUdim1_Ud0 + (1+a)*dt^2*(.5-b)*dUddim1_Ud0) +...
            dFnldud*(dUddim1_Ud0 + (1+a)*dt*(1-g)*dUddim1_Ud0) -...
            dFnldu0*dUim1_Ud0 - dFnldud0*dUdim1_Ud0);
        
        dUddi_U0 = -J\dRdU0;
        dUddi_Ud0 = -J\dRdUd0;
        
        dUi_U0 = dUim1_U0 + (1+a)*dt*dUdim1_U0 + (1+a)*dt^2*((.5-b)*dUddim1_U0 + b*dUddi_U0);
        dUi_Ud0 = dUim1_Ud0 + (1+a)*dt*dUdim1_Ud0 + (1+a)*dt^2*((.5-b)*dUddim1_Ud0 + b*dUddi_Ud0);
        
        dUdi_U0 = dUdim1_U0 + (1+a)*dt*((1-g)*dUddim1_U0 + g*dUddi_U0);
        dUdi_Ud0 = dUdim1_Ud0 + (1+a)*dt*((1-g)*dUddim1_Ud0 + g*dUddi_Ud0);

	    % Propagate
        dUim1_U0 = dUi_U0;
        dUim1_Ud0 = dUi_Ud0;
        
        dUdim1_U0 = dUi_U0;
        dUdim1_Ud0 = dUi_Ud0;
        
        dUddim1_U0 = dUddi_U0;
        dUddim1_Ud0 = dUddi_Ud0;
        
        dFnldu0 = dFnldu;
        dFnldud0 = dFnldud;
      end
      
      if strcmp(opts.Display, 'progress') || strcmp(opts.Display, 'both')
          fprintf('%d: %.4e/%.4e %.4e\n', i, T(i), T1, dt);
          if strcmp(opts.Display, 'both')
            fprintf('---------------------------------------------------\n');
          end
      end
      
      %% Check for kill
      if strcmp(opts.Display, 'waitbar')
          waitbar(T(i)/T1, wb, sprintf('Progress: %.4e/%.4e dt: %.4e', T(i), T1, dt))
          
          if (~ishandle(wb))
            break;
          elseif getappdata(wb, 'interrupt')
            delete(wb);
       
            U   =   U(:, 1:i);
            Ud  =  Ud(:, 1:i);
            Udd = Udd(:, 1:i);
            T   =   T(1:i);
            
            if nargout==11
              varargout{1} = dUi_U0;
              varargout{2} = dUi_Ud0;
              varargout{3} = dUdi_U0;
              varargout{4} = dUdi_Ud0;
              varargout{5} = dUddi_U0;
              varargout{6} = dUddi_Ud0;
            end
            
            return;
          end
      end
  end
  if strcmp(opts.Display, 'waitbar')  
    waitbar(1.0, wb, 'COMPLETED!');
    delete(wb);
  end
  
  if nargout==11
    varargout{1} = dUi_U0;
    varargout{2} = dUi_Ud0;
    varargout{3} = dUdi_U0;
    varargout{4} = dUdi_Ud0;
    varargout{5} = dUddi_U0;
    varargout{6} = dUddi_Ud0;
  end
end

% Command to delete stray waitbars: delete(findall(0,'type','figure','tag','TMWWaitbar'))