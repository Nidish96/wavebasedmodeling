function [Me, Ke] = PLANTMBMEMATS(Le, pars)
  
  rho = pars.rho;
  E = pars.E;
  G = pars.G;
  Ar = pars.A;
  Iz = pars.Iz;
  
%  % 1D bar (axial)
%  Me_bar = rho*Ar*Le/6*[2 1;
%                       1 2];
%  Ke_bar = E*Ar/Le*[1 -1;
%                   -1 1];
  
%  % Planar timoshenko (transverse)
%  Me_tb = rho*Le/6*[2*Ar 0 Ar 0;
%	                  0 2*Iz 0 Iz;
%	                  Ar 0 2*Ar 0;
%	                  0 Iz 0 2*Iz];

%  Ke_tb = 1/Le*[G*Ar 0 -G*Ar 0;
%	              0 E*Iz 0 -E*Iz;
%	             -G*Ar 0 G*Ar 0;
%	              0 -E*Iz 0 E*Iz] +...
%          G*Ar*Le/4*[0 0 0 0;
%		                 0 1 0 1;
%		                 0 0 0 0;
%		                 0 1 0 1] - ...
%          G*Ar/2*[0 1 0 1;
%	                1 0 -1 0;
%	                0 -1 0 -1;
%	                1 0 -1 0];
  
%  Me = zeros(6);
%  Ke = zeros(6);
  
%  Me([1 4], [1 4]) = Me_bar;
%  Me([2 3 5 6], [2 3 5 6]) = Me_tb;
  
%  Ke([1 4], [1 4]) = Ke_bar;
%  Ke([2 3 5 6], [2 3 5 6]) = Ke_tb;
  
  Me = rho*Le/6*[2*Ar, 0   , 0   , Ar  , 0   , 0;
                 0   , 2*Ar, 0   , 0   , Ar  , 0;
                 0   , 0   , 2*Iz, 0   , 0   , Iz;
                 Ar  , 0   , 0   , 2*Ar, 0   , 0;
                 0   , Ar  , 0   , 0   , 2*Ar, 0;
                 0   , 0   , Iz  , 0   , 0   , 2*Iz];
  Ke = 1/Le*[E*Ar , 0    , 0    , -E*Ar, 0    , 0;
             0    , G*Ar , 0    , 0    , -G*Ar, 0;
             0    , 0    , E*Iz , 0    , 0    , -E*Iz;
             -E*Ar, 0    , 0    , E*Ar , 0    , 0;
             0    , -G*Ar, 0    , 0    , G*Ar , 0;
             0    , 0    , -E*Iz, 0    , 0    , E*Iz] + ...
       G*Ar*Le/4*[0 0 0 0 0 0;
                  0 0 0 0 0 0
                  0 0 1 0 0 1;
                  0 0 0 0 0 0;
                  0 0 0 0 0 0;
                  0 0 1 0 0 1] + ...
       G*Ar/2*[0 , 0 , 0 , 0 , 0 , 0;
               0 , 0 , 1 , 0 , 0 , 1;
               0 , 1 , 0 , 0 , -1, 0;
               0 , 0 , 0 , 0 , 0 , 0;
               0 , 0 , -1, 0 , 0 , -1;
               0 , 1 , 0 , 0 , -1, 0];
end
