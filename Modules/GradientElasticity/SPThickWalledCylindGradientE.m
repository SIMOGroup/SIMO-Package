function ur = SPThickWalledCylindGradientE(ri,ro,pi,po,nuy,l)
% S = SPThickWalledCylindExactStress(Ri, Ro, pi, po, nu)
% Evaluate analytical solution of stresses for "Thick-Walled Cylinder"
% problem
% ------------------------------------------------------------------
% Input:
%       ri: inner radius
%       ro: outer radius
%       pi: inner pressure
%       po: outer pressure
%       nuy: Poisson's ratio
%       l: internal length scale
% ------------------------------------------------------------------
% Output:
%       ur: dispalcement in R direction
% ------------------------------------------------------------------

%{
Copyright (C) <2014-2016>  <Thai Son, Khanh Chau-Nguyen, Hung Nguyen-Xuan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

lamda = 7000;
muyt = 3000;
a1=0;
a2=0.5*l^2*lamda;
a3=0;
a4=1^2*muyt;
a5=0;

l_head = sqrt(2*(a1+a2+a3+a4+a5)/(lamda+2*muyt));

xi = (a1+a2+a3+a4+a5)/(2*(a4+a5));
zet = (lamda+muyt)/muyt;

phi_a = (ri^3*xi+ri)*besseli(1,ri)-(ri^2)*besselj(1,ri);
psi_a = (ri^3*xi+ri)*besselk(1,ri)-(ri^2)*bessely(1,ri);

phi_b = (ro^3*xi+ro)*besseli(1,ro)-(ri^2)*besselj(1,ri);
psi_b = (ro^3*xi+ro)*besselk(1,ro)-(ro^2)*besseli(1,ri);

chi = -(1/zet)*(1/ro^2-1/ri^2)*(phi_a*psi_b-phi_b*psi_a)-(besseli(1,ro)/ro...
    - besseli(1,ri)/ri)*(psi_b-psi_a)+(besselk(1,ro)/ro-besselk(1,ri)/ri)*(phi_b-phi_a);

C1 = po/(2*(lamda+muyt)*chi)*(1/(zet*ri^2)*(phi_a*psi_b-phi_b*psi_a)+besseli(1,ri)/ri*(psi_b-psi_a)-...
    besselk(1,ri)/ri*(phi_b-phi_a));
C2 = po/(2*(lamda+muyt))*(phi_a*psi_b-phi_b*psi_a)/chi;
C3 = -po/(2*(lamda+muyt))*(psi_b-psi_a)/chi;
C4 = -po/(2*(lamda+muyt))*(phi_b-phi_a)/chi;

%r = linspace(ri, ro, 11);
%rho = (r/l_head);
ur.disp = @(r) l_head*(C1.*(r/l_head)+C2./(r/l_head)+C3*besseli(1,(r/l_head))+C4*besselk(1,(r/l_head)));
end



