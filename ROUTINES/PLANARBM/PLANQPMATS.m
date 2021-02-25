function [Q, T, Xqps] = PLANQPMATS(Nq, ypos, bm, pars)
%PLANQPMATS returns the quadrture oisition and integration matrices for planar timoshenko beams with DOFs assumed to be [ux; uy; thetaz].
%
%  USAGE:
%    [Q, T, Xqps] = PLANQPMATS(Nq, ypos, bm, pars);
%    % Q*unodal gives the [ux;uy;thz] at quadrature locations
%    % T'*tqp integrates the quadrature point tractions to nodal forces
%  INPUTS:
%    Nq		: Quadrature points per element
%    ypos	: Location on the plane (relative to neutral axis) of quadrature points 
%    bm		: Beam structure with, fields
%	X, Y	: (Ne+1,1) vectors of X and Y poisitions of the nodes
%  	Wi,Wo   : In-plane and out-of-plane widths
%  OUTPUTS:
%    Q, T 	: ((Ne*Nq)*3, (Ne+1)*3) matrices as above
%    Xqps	: (Ne*Nq,1) X-locations of the quadrature points
  Ne = length(bm.X)-1;

  Q = zeros((Ne*Nq)*3, (Ne+1)*3);
  T = zeros((Ne*Nq)*3, (Ne+1)*3);
  Xqps = zeros(Ne*Nq, 1);
  for e = 1:Ne
    n1 = e;
    n2 = e+1;

    is = (n1-1)*3+1;
    ie = n2*3;

    Le = bm.X(n2)-bm.X(n1);

    [xi, wi] = LGWT(Nq, 0, Le);
    xi = xi(end:-1:1);
    wi = wi(end:-1:1);
    Xqps((e-1)*Nq+(1:Nq)) = bm.X(n1) + xi;

    wi = wi*bm.Wo;  % Out of plane width
    for qi = 1:Nq
      Nsf = [1-xi(qi)/Le, 0, 0, xi(qi)/Le, 0, 0;
	     0, 1-xi(qi)/Le, 0, 0, xi(qi)/Le, 0;
	     0, 0, 1/2, 0, 0, 1/2];
      
      Q((e-1)*(Nq*3)+(qi-1)*3+(1:3), is:ie) = [1 0 -ypos;0 1 0;0 0 1]*Nsf;

      T((e-1)*(Nq*3)+(qi-1)*3+(1:3), is:ie) = [1 0 -ypos;0 1 0;0 0 1]*Nsf*wi(qi);
    end
  end
end
