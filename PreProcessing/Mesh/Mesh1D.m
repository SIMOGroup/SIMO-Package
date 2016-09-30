function Mesh = Mesh1D(NURBS)
% function Mesh = Mesh1D(NURBS)
% Create mesh structure for 1D problems
% -------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
% -------------------------------------------------------------------
% Output:
%       Mesh: Mesh structure in struct format, where
%           Mesh.Dof -- degree of freedom per node
%           Mesh.NEN -- number of local basis functions
%           Mesh.NEl -- number of elements
%           Mesh.El -- element connectivity
%           Mesh.NDof -- number of degree of freedom
%           Mesh.Boundary(1).Dofs = dofs index at side 1 (u = 0)
%           Mesh.Boundary(2).Dofs = dofs index at side 2 (u = 1)
% -------------------------------------------------------------------

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

El = BuildConn(NURBS.Order, NURBS.KntVect{1});

Mesh.Dof = 1;
Mesh.NEN = NURBS.Order + 1; % number of local basis functions
Mesh.NEl = numel(unique(NURBS.KntVect{1})) - 1;
Mesh.El = El;
Mesh.NDof = NURBS.NNP;
Mesh.Boundary(1).Dofs = 1;
Mesh.Boundary(2).Dofs = NURBS.NCtrlPts(1);
end
