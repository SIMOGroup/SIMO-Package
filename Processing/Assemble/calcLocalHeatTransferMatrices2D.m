function [HVals, FVals] = calcLocalHeatTransferMatrices2D(NURBS, Mesh, h, Ta, Ref)
% function [HVals, FVals] = calcLocalHeatTransferMatrices2D(NURBS, Mesh, h, Ta, Ref)
% -----------------------------------------------------------------
% Evaluate local surface heat transfer matrices for 2D
% -----------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       Mesh: Mesh structure
%       h: convective heat transfer coefficient
%       Ta: ambient temperature
%       Ref: referenced index to indicate boundary curve,
%       ******************************************************
%       *                                                    *
%       *                   Refs(i) = 4                      *
%       *                   (v = 1)                          *
%       *               o------------o              v        *
%       *               |            |              ^        *
%       *               |            |              |        *
%       *   Refs(i) = 1 |            | Refs(i) = 2  |        *
%       *       (u = 0) |            |   (u = 1)    +----> u *
%       *               |            |                       *
%       *               o------------o                       *
%       *                  Refs(i) = 3                       *
%       *                   (v = 0)                          *
%       ******************************************************
% --------------------------------------------------------------------
% Output:
%       HVals: matrix stores surface heat transfer matrices
%           (size(HVals) = [NEN ^ 2, NEL])
%       FVals: matrix stores local body forces
%           (size(FVals) = [NEN, NEL])
% --------------------------------------------------------------------

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

MeshBdry = Mesh.Boundary(Ref);
NURBSBdry = NURBSBoundary(NURBS, Ref);
NGPs = NURBSBdry.Order + 1;
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBSBdry.Order, NURBSBdry.NCtrlPts, NURBSBdry.KntVect{1}, 1, NGPs, MeshBdry.NElDir);
CtrlPts = NURBSBdry.CtrlPts3D(1 : 2, :)';
Weights = NURBSBdry.Weights;

HVals = zeros(MeshBdry.NEN ^ 2, MeshBdry.NEl);
FVals = zeros(MeshBdry.NEN, MeshBdry.NEl);
for ex = 1 : MeshBdry.NElDir
    HL = zeros(MeshBdry.NEN); % local surface heat transfer matrix
    FL = zeros(MeshBdry.NEN, 1);
    for qx = 1 : NGPs
        
        N0 = reshape(Nx(ex, qx, :, 1), 1, []);
        N1 = reshape(Nx(ex, qx, :, 2), 1, []);
        
        [R0, R1] = Rationalize(Weights(MeshBdry.El(ex, :)), N0, N1);
        
        dxdxi = R1 * CtrlPts(MeshBdry.El(ex, :), :);
        
        J1 = norm(dxdxi);
        
        HL = HL + R0' * h * R0 * J1 * Jx(ex) * Wx(qx);
        
        FL = FL + R0' * h * Ta * J1 * Jx(ex) * Wx(qx);
    end
    FVals(:, ex) = FL;
    HVals(:, ex) = HL(:);
end
end