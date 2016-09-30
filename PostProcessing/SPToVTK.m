function SPToVTK(NURBS, d, ParaPts, filename, fieldname)
% SPToVTK(NURBS, d, ParaPts, filename, fieldname)
% --------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       d: temperature or displacement field
%       ParaPts: parameter points per direction
%       filename: name of file
%       fieldname: name of field
% ------------------------------------------------------------------
% Output:
%       filename.vts
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

F1 = reshape(d, NURBS.NNP, [])';
F2 = reshape(F1, [size(F1, 1), NURBS.NCtrlPts]);
F2w = bsxfun(@times, F2, NURBS.Weights);
CFw = cat(1, NURBS.CtrlPts4D, F2w);
tmp = BsplineEval(NURBS.KntVect, CFw, ParaPts);
Weights = tmp(4, :, :, :);
C = bsxfun(@rdivide, tmp(1 : 3, :, :, :), Weights);
F = bsxfun(@rdivide, tmp(5 : end, :, :, :), Weights);
exportToVTK(C, squeeze(F), filename, fieldname);
end