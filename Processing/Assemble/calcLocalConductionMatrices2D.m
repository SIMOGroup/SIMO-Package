function [KVals, FVals] = calcLocalConductionMatrices2D(Mesh, NURBS, ka, s)
% [KVals, FVals] = calcLocalConductionMatrices2D(Mesh, NURBS, ka, s)
% --------------------------------------------------------------------
% Evaluate local conduction matrices and body forces
% --------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       ka: conductivity coefficient
%       s: heat source
% --------------------------------------------------------------------
% Output:
%       KVals: matrix stores local conduction matrices
%           (size(KVals) = [NEN ^ 2, NEL])
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

% conduction matrix
D = ka * eye(2); % D = k * ones(NSD, NSD), NSD: number of space dimension

NGPs = NURBS.Order + 1;

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));

N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);

KVals = zeros(Mesh.NEN ^ 2, Mesh.NEl);
FVals = zeros(Mesh.NEN, Mesh.NEl);

Weights = reshape(NURBS.Weights, 1, []);
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';
for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        e = sub2ind(Mesh.NElDir, ex, ey);
        Ke = zeros(Mesh.NEN);
        Fe = zeros(Mesh.NEN, 1);
        for qy = 1 : NGPs(2)
            for qx = 1 : NGPs(1)
                k = 1;
                for j = 1 : NURBS.Order(2) + 1
                    for i = 1 : NURBS.Order(1) + 1
                        N0(k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1);
                        N1(1, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 1);
                        N1(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 2);
                        k = k + 1;
                    end
                end
                [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0, N1);
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % Gradient of mapping from parameter space to physical 
                % space
                dxdxi = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
                
                J1 = det(dxdxi);
                
                % Compute derivatives of basis functions w.r.t
                % physical coordinates
                dRdx = dxdxi^(-1) * R1;
                
                B = dRdx; % size(B) = [NSD, NEN]
                
                % compute element stiffness at quadrature point
                Ke = Ke + B' * D * B * J1 * J2 * W;
                
                Pts = R0 * CtrlPts(Mesh.El(e, :), :);
                
                x = Pts(1); y = Pts(2);
                
                % element source vector
                Fe = Fe + R0' * s(x, y) * J1 * J2 * W; 
            end
        end
        FVals(:, e) = Fe;
        KVals(:, e) = Ke(:);
    end
end
end