function Cw = NURBSLine(P1, P2)
% Cw = NURBSLine(P1, P2)
% -----------------------------------------------------------------
% Contruct a straight line
% -----------------------------------------------------------------
% Input:
%       P1: the 1st point of the line
%       P2: the 2nd point of the line
% -----------------------------------------------------------------
% Output:
%       Cw: NURBS structure of the line
% -----------------------------------------------------------------

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

CtrlPts = zeros(4, 2);
CtrlPts(1 : numel(P1), 1) = P1;
CtrlPts(1 : numel(P2), 2) = P2;
CtrlPts(4, :) = 1;
KntVect = [0 0 1 1];
Cw = CreateNURBS({KntVect}, CtrlPts);
end