function ONURBS = NURBSRepara(INURBS, Dir, StartKnt, EndKnt)
% function ONURBS = NURBSRepara(INURBS, Dir, StartKnt, EndKnt)
% ------------------------------------------------------------------
% Reparameterize the knot vectors
% ------------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       Dir: direction to reparameter
%           (xi--Dir = 1, eta--Dir = 2, zeta--Dir = 3)
%       StartKnt: starting knot value
%       EndKnt: ending knot value
% ------------------------------------------------------------------
% Output:
%       ONURBS: output NURBS
% ------------------------------------------------------------------

%{
Copyright (C) <2014-2016>  <Khanh Chau-Nguyen, Hung Nguyen-Xuan>

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

p = INURBS.Order(Dir);
KntVect = INURBS.KntVect{Dir};
a = KntVect(p + 1);
b = KntVect(end - p);
if (StartKnt == a && EndKnt == b)
    ONURBS = INURBS;
else
    assert(StartKnt < EndKnt, 'StartKnot must be less than EndKnot')
    
    KntVect = KntVect - a;
    KntVect = KntVect * (EndKnt - StartKnt) / (b - a);
    KntVect = KntVect + StartKnt;
    
    finKntVect = INURBS.KntVect;
    finKntVect{Dir} = KntVect;
    
    ONURBS = CreateNURBS(finKntVect, INURBS.CtrlPts4D);
end
end