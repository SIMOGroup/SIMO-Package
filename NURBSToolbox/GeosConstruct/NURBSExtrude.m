function ONURBS = NURBSExtrude(INURBS, DisplVect)
% ONURBS = NURBSExtrude(INURBS, DisplVect)
% --------------------------------------------------------------------
% Construct a NURBS surface (volume) by extruding a NURBS curve (surface)
% ----------------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       DisplVect: displacment vector
%       (e.g. DisplVect = [0 0 11] means the height of the extrusion
%              is 11 units along the z axis)
% ----------------------------------------------------------------------
% Ouput:
%       ONURBS: input NURBS structure
% ----------------------------------------------------------------------

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

% translation matrix
A = eye(4);
A(1 : end - 1, 4) = DisplVect';

a = INURBS.CtrlPts4D;
b = reshape(A * reshape(a, 4, []), size(a));

CtrlPts = cat(INURBS.Dim + 2, a, b);

ONURBS = CreateNURBS([INURBS.KntVect, [0, 0, 1, 1]], CtrlPts);
end