function srf = NURBSRuled(crv1, crv2)

% NRBRULED: Construct a ruled surface between two NURBS curves.
% 
% Calling Sequence:
% 
%   srf = nrbruled(crv1, crv2)
% 
% INPUT:
% 
%   crv1	: First NURBS curve, see nrbmak.
% 
%   crv2	: Second NURBS curve, see nrbmak.
%
% OUTPUT:
% 
%   srf		: Ruled NURBS surface.
% 
% Description:
% 
%   Constructs a ruled surface between two NURBS curves. The ruled surface is
%   ruled along the V direction.
% 
% Examples:
% 
%   Construct a ruled surface between a semicircle and a straight line.
% 
%   cir = nrbcirc(1,[0 0 0],0,pi);
%   line = nrbline([-1 0.5 1],[1 0.5 1]);
%   srf = nrbruled(cir,line);
%   nrbplot(srf,[20 20]);
%
%    Copyright (C) 2000 Mark Spink
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% ensure both curves have a common degree
d = max ([crv1.Order, crv2.Order]);
crv1 = PRefine(crv1, 1, d - crv1.Order);
crv2 = PRefine(crv2, 1, d - crv2.Order);

% merge the knot vectors, to obtain a common knot vector
k1 = crv1.KntVect{1};
k2 = crv2.KntVect{1};
ku = unique ([k1 k2]);
n = length (ku);
ka = [];
kb = [];
for i = 1:n
  i1 = length (find (k1 == ku(i)));
  i2 = length (find (k2 == ku(i)));
  m = max (i1, i2);
  ka = [ka ku(i)*ones(1,m-i1)];
  kb = [kb ku(i)*ones(1,m-i2)];
end
crv1 = HRefine(crv1, 1, ka);
crv2 = HRefine(crv2, 1, kb);

CtrlPts(:, :, 1) = crv1.CtrlPts4D;
CtrlPts(:, :, 2) = crv2.CtrlPts4D;
srf = CreateNURBS({crv1.KntVect{1}, [0 0 1 1]}, CtrlPts);
end