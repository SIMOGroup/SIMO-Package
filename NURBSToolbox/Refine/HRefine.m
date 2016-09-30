function ONURBS = HRefine(INURBS, Dir, InsKnts)
% ONURBS = HRefine(INURBS, Dir, InsKnts)
% --------------------------------------------------------------
% Knot refinement for curve, surface and solid
%---------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       Dir: parameter direction to refine, xi = 1, eta = 1...
%       InsKnts: a list of knots to insert
%---------------------------------------------------------------
% Output:
%       ONURBS: output NURBS structure
% --------------------------------------------------------------

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

if isempty(InsKnts)
    ONURBS = INURBS;
else
    % if 'Dir'th dimension of CtrlPts4D matrix is not the last
    % dimension, then permute it until it lies at last dimension
    if Dir ~= INURBS.Dim
        Dirs = 1 : INURBS.Dim + 1;
        Dirs(Dir + 1) = [];
        Dirs = [Dirs, Dir + 1];
        tmp = permute(INURBS.CtrlPts4D, Dirs);
    else
        tmp = INURBS.CtrlPts4D;
    end
    dim = size(tmp);
    tmp = reshape(tmp, [], dim(end));
    
    [KntVect, CtrlPts] = RefineKntVectCurv(INURBS.Order(Dir), INURBS.KntVect{Dir}, tmp, InsKnts);
    
    CtrlPts = reshape(CtrlPts, [dim(1 : end - 1), size(CtrlPts, 2)]);
    
    if Dir ~= INURBS.Dim
        % permute the last dimension to initial position. The 
        % positions of the other dimensions do not change 
        % relative to one another.
        Dirs = 1 : INURBS.Dim + 1;
        DirsTmp = Dirs(Dir + 1 : end - 1);
        Dirs(Dir + 1) = INURBS.Dim + 1;
        Dirs(Dir + 2 : end) = DirsTmp;        
        CtrlPts = permute(CtrlPts, Dirs);
    end
    
    OKntVect = INURBS.KntVect;
    OKntVect{Dir} = KntVect;
    ONURBS = CreateNURBS(OKntVect, CtrlPts);
end
end
