function ONURBS = NURBSJoin(INURBS1, INURBS2, Dir)
% function ONURBS = NURBSJoin(INURBS1, INURBS2, Dir)
% -------------------------------------------------------------
% Joint two NURBS structures
% Notice: only use to join two coarsest NURBS structures,
% not for refined NURBS
% -------------------------------------------------------------
% Input:
%       INURBS1: first NURBS structure
%       INURBS2: second NURBS structure
%       Dir: direction to join (xi: Dir = 1, eta: Dir = 2, 
%       zeta: Dir = 3)
% -------------------------------------------------------------
% Output:
%       ONURBS: ouput NURBS structure
% -------------------------------------------------------------

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

dim = INURBS1.Dim;
assert(dim == INURBS2.Dim);
assert(Dir >=1 && Dir <= dim)

order(1) = INURBS1.Order(Dir);
order(2) = INURBS2.Order(Dir);

IKnts{1} = INURBS1.KntVect{Dir};
IKnts{2} = INURBS2.KntVect{Dir};

Knt(1) = IKnts{1}(end - order(1));
Knt(2) = IKnts{2}(order(2) + 1);
Knt(3) = IKnts{2}(end - order(2));

INURBS2 = NURBSRepara(INURBS2, Dir, Knt(2) + Knt(1), Knt(3) + Knt(1));

p = max(order);
INURBS1 = PRefine(INURBS1, Dir, p - order(1));
INURBS2 = PRefine(INURBS2, Dir, p - order(2));

CtrlPts{1} = INURBS1.CtrlPts4D;
CtrlPts{2} = INURBS2.CtrlPts4D;

n1 = size(CtrlPts{1});
idx1 = cell(1, numel(n1));
for i = 1 : numel(n1)
   idx1{i} = 1 : n1(i); 
end
idx1{Dir + 1} = 1 : n1(Dir + 1) - 1;

n2 = size(CtrlPts{2});
idx2 = cell(1, numel(n2));
for i = 1 : numel(n2)
   idx2{i} = 1 : n2(i); 
end
idx2{Dir + 1} = 2 : n2(Dir + 1);

ExCtrlPts{1} = CtrlPts{1}(idx1{:}); % extracted control points
ExCtrlPts{3} = CtrlPts{2}(idx2{:});

idx1{Dir + 1} = n1(Dir + 1);
idx2{Dir + 1} = 1;

ExCtrlPts{2} = (CtrlPts{1}(idx1{:}) + CtrlPts{2}(idx2{:})) / 2;

A = cat(Dir + 1, ExCtrlPts{:});

IKnts{1} = INURBS1.KntVect{Dir};
IKnts{2} = INURBS2.KntVect{Dir};

OKnts{1} = IKnts{1}(1 : end - p - 1);
OKnts{2} = ones(1, p) * Knt(1);
OKnts{3} = IKnts{2}(p + 2 : end);

OKnts = cell2mat(OKnts);
OKnts = OKnts ./ max(OKnts);

KntVect = INURBS1.KntVect; % knot vectors
KntVect{Dir} = OKnts;

ONURBS = CreateNURBS(KntVect, A);

end