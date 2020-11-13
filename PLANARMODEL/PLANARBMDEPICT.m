function [] = PLANARBMDEPICT(U, BM, varargin)
  
  col = 'r';
  alph = 0.1;
  if nargin>=3
    col = varargin{1};
  end 
  if nargin>=4
    alph = varargin{2};
  end
  
%% Small Angle Assumption  
%   Top_X = BM.X+U(1:3:end) - U(3:3:end)*BM.Wi/2;
%   Top_Y = BM.Y+U(2:3:end) + BM.Wi/2;
%   
%   Bot_X = BM.X+U(1:3:end) + U(3:3:end)*BM.Wi/2;
%   Bot_Y = BM.Y+U(2:3:end) - BM.Wi/2;  
  
%% Large Angle  
 Top_X = BM.X+U(1:3:end) - sin(U(3:3:end))*BM.Wi/2;
 Top_Y = BM.Y+U(2:3:end) + cos(U(3:3:end))*BM.Wi/2;
  
 Bot_X = BM.X+U(1:3:end) + sin(U(3:3:end))*BM.Wi/2;
 Bot_Y = BM.Y+U(2:3:end) - cos(U(3:3:end))*BM.Wi/2;
  
  
  %% Plotting
  plot(Top_X, Top_Y, 'k.-'); hold on 
  plot(Bot_X, Bot_Y, 'k.-'); hold on 
  
  for e=1:(length(BM.X)-1)
    X = [Top_X([e e+1]); Bot_X([e+1 e])];
    Y = [Top_Y([e e+1]); Bot_Y([e+1 e])];
    
    n1 = e;
    n2 = e+1;
    
    is = (n1-1)*3+1;
    ie = n2*3;
    
    Le = BM.X(e+1)-BM.X(e);
    ex = [ones(4,1) zeros(4,1) -Y]*kron([-1 1]/Le, eye(3))*U(is:ie);
    gxy = [[zeros(4,1) ones(4,1) zeros(4,1)]*kron([-1 1]/Le, eye(3)) + [zeros(4,2) -ones(4,1)]*kron([1 1]/2, eye(3))]*U(is:ie);
    
    p1 = (ex+sqrt(ex.^2+4*gxy.^2))/2;
    p2 = (ex-sqrt(ex.^2+4*gxy.^2))/2;

    if strcmp(col, 'ex')
      fill(X, Y, ex, 'FaceAlpha', alph); hold on 
    elseif strcmp(col, 'gxy')
      fill(X, Y, gxy, 'FaceAlpha', alph); hold on 
    elseif strcmp(col, 'p1')
      fill(X, Y, p1, 'FaceAlpha', alph); hold on 
    elseif strcmp(col, 'p2')
      fill(X, Y, p2, 'FaceAlpha', alph); hold on 
    else 
      fill(X, Y, col, 'FaceAlpha', alph); hold on 
    end
  end
  if strcmp(col, 'ex') || strcmp(col, 'gxy')|| strcmp(col, 'p1')|| strcmp(col, 'p2')
    plot(BM.X+U(1:3:end), BM.Y+U(2:3:end), 'ko-'); hold on 
  else 
    plot(BM.X+U(1:3:end), BM.Y+U(2:3:end), 'o-', 'Color', col); hold on 
  end 
end
