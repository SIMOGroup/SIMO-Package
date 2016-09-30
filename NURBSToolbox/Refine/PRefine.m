function ONURBS = PRefine(INURBS, Dir, t)
% function ONURBS = PRefine(INURBS, Dir, t)
% --------------------------------------------------------------
% Order elevation for curve, surface and solid
%---------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       Dir: parameter direction to refine, xi = 1, eta = 2...
%       t: times
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

if (t == 0)
    ONURBS = INURBS;
else
    % if 'Dir'th dimension of CtrlPts matrix is not the last
    % dimension, then permute it until it lies at the last dimension
    if Dir ~= INURBS.Dim
        Dirs = 1 : INURBS.Dim + 1;
        Dirs(Dir + 1) = [];
        Dirs = [Dirs, Dir + 1];
        temp = permute(INURBS.CtrlPts4D, Dirs);
    else
        temp = INURBS.CtrlPts4D;
    end
    dim = size(temp);
    temp = reshape(temp, [], dim(end));
    
    [KntVect, CtrlPts] = DegreeElevateCurv(...
        INURBS.Order(Dir), INURBS.KntVect{Dir}, temp, t);
    
    CtrlPts = reshape(CtrlPts, [dim(1 : end - 1), size(CtrlPts, 2)]);
    
    if Dir ~= INURBS.Dim
        % permute the last dimension to initial position. The 
        % positions of the other dimensions do not change 
        % relative to one another.
        Dirs = 1 : INURBS.Dim + 1;
        DirsTemp = Dirs(Dir + 1 : end - 1);
        Dirs(Dir + 1) = INURBS.Dim + 1;
        Dirs(Dir + 2 : end) = DirsTemp;        
        CtrlPts = permute(CtrlPts, Dirs);
    end
    
    finKntVect = INURBS.KntVect;
    finKntVect{Dir} = KntVect;
    ONURBS = CreateNURBS(finKntVect, CtrlPts);
end
end