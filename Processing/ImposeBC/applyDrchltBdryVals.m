function [f, d, FreeIdcs] = applyDrchltBdryVals(BdryIdcs, BdryVals, K, f)
% function [f, d, freeIdcs] = applyDrchltBdryVals(BdryIdcs, BdryVals, K, f)
% Apply boundary values for Dirichlet boundary condition
% ----------------------------------------------------------------
% Input: 
%       BdryIdcs: indices of boundary control points
%       BdryVals: corresponding values at these control points
%       K: global stiffness matrix
%       f: force vector (right hand side)
% Output:
%       f: modified force vector
%       d: displacement vector
%       freeIdcs: vector contains indices of unknowns
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

NDof = numel(f);
FreeIdcs = setdiff(1 : NDof, BdryIdcs);
d = zeros(NDof, 1);
d(BdryIdcs) = BdryVals;
f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 
end