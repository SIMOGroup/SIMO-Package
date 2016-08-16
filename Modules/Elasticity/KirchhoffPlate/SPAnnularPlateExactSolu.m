function solu = SPAnnularPlateExactSolu(R1, R0, t, E, nu, q0)
% function solu = SPAnnularPlateExactSolu(R1, R0, t, E, nu, q0)
% Input:
%       R1: inner radius
%       R0: outer radius
%       t: thickness
%       E: Young's modulus
%       nu: Poisson's ratio
%       q0: uniformly distributed force
% Outut:
%       solu.w: displacement
%       solu.theta: rotation
% ----------------------------------------------------------------------

r = @(x, y) sqrt(x.^2 + y.^2);
D = E*t^3/(12*(1-nu^2));

C1 = -((1/32)*R0^2-(1/32)*R1^2)*q0*(4*log(R0)*(R0^2)*(R1^2)-4*log(R1)*(R0^2)*(R1^2)-R0^4+R1^4)*(1/D)*(1/(4*log(R0)^2*(R0^2)*(R1^2)-8*log(R0)*log(R1)*(R0^2)*(R1^2)+4*log(R1)^2*(R0^2)*(R1^2)-R0^4+2*R0^2*(R1^2)-R1^4));
C2 = -(1/64)*q0*(8*log(R1)^2*(R0^4)*(R1^2)-8*log(R1)*log(R0)*(R0^4)*(R1^2)-8*log(R1)*log(R0)*(R0^2)*(R1^4)+8*log(R0)^2*(R0^2)*(R1^4)+2*log(R1)*(R0^4)*(R1^2)-4*log(R1)*(R0^2)*(R1^4)+2*log(R1)*(R1^6)+2*log(R0)*(R0^6)-4*log(R0)*(R0^4)*(R1^2)+2*log(R0)*(R0^2)*(R1^4)-R0^6+R0^4*(R1^2)+R0^2*(R1^4)-R1^6)*(1/D)*(1/(4*log(R1)^2*(R0^2)*(R1^2)-8*log(R1)*log(R0)*(R0^2)*(R1^2)+4*log(R0)^2*(R0^2)*(R1^2)-R0^4+2*R0^2*(R1^2)-R1^4));
C3 = -(1/16)*R0^2*q0*(R1^2)*(R0^4*log(R1)-log(R1)*(R1^4)-R0^4*log(R0)+log(R0)*(R1^4)+R0^4-2*R0^2*(R1^2)+R1^4)*(1/D)*(1/(4*log(R1)^2*(R0^2)*(R1^2)-8*log(R1)*log(R0)*(R0^2)*(R1^2)+4*log(R0)^2*(R0^2)*(R1^2)-R0^4+2*R0^2*(R1^2)-R1^4));
C4 = (1/64)*R0^2*q0*(R1^2)*(4*log(R1)^2*(R0^4)-4*log(R1)*log(R0)*(R0^4)-4*log(R1)*log(R0)*(R1^4)+4*log(R0)^2*(R1^4)+2*R0^4*log(R1)-4*log(R1)*(R0^2)*(R1^2)+2*log(R1)*(R1^4)+2*R0^4*log(R0)-4*log(R0)*(R0^2)*(R1^2)+2*log(R0)*(R1^4)-R0^4+2*R0^2*(R1^2)-R1^4)*(1/D)*(1/(4*log(R1)^2*(R0^2)*(R1^2)-8*log(R1)*log(R0)*(R0^2)*(R1^2)+4*log(R0)^2*(R0^2)*(R1^2)-R0^4+2*R0^2*(R1^2)-R1^4));

solu.w = @(x, y) C1*r(x,y).^2.*log(r(x,y))+C2*r(x,y).^2+C3*log(r(x,y))+C4+(q0*r(x,y).^4)/(64*D);
solu.theta = @(x, y) 2*C1*r(x,y).*log(r(x,y))+C1*r(x,y)+2*C2*r(x,y)+C3./r(x,y)+(1/16)*q0*r(x,y).^3/D;

end
