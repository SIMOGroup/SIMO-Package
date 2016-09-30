function [KVals, FVals] = calcLocalConductionMatrices3D(Mesh, NURBS, kappa, s)
% function [KVals, FVals] = calcLocalConductionMatrices3D(Mesh, NURBS, kappa, s)
% -----------------------------------------------------------------------
% Evaluate local conduction matrices and body forces
% --------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       kappa: conductivity coefficient
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

% conductivity matrix
D = kappa * eye(3); % D = k * ones(NSD, NSD), NSD: number of space dimension

NGPs = NURBS.Order + 1;

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));
[Jz, Wz, ~, Nz] = calcDersBasisFunsAtGPs(NURBS.Order(3), NURBS.NCtrlPts(3), NURBS.KntVect{3}, 1, NGPs(3), Mesh.NElDir(3));

N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);

KVals = zeros(Mesh.NEN ^ 2, Mesh.NEl);
FVals = zeros(Mesh.NEN, Mesh.NEl);

CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';
Weights = reshape(NURBS.Weights, 1, []);

for ez = 1 : Mesh.NElDir(3)
    for ey = 1 : Mesh.NElDir(2)
        for ex = 1 : Mesh.NElDir(1)
            e = sub2ind(Mesh.NElDir, ex, ey, ez);
            ElConn = Mesh.El(e, :); %  element connectivity
            Ke = zeros(Mesh.NEN, Mesh.NEN);
            Fe = zeros(Mesh.NEN, 1);
            for qz = 1 : NGPs(3)
                for qy = 1 : NGPs(2)
                    for qx = 1 : NGPs(1)
                        l = 1;
                        for k = 1 : NURBS.Order(3) + 1
                            for j = 1 : NURBS.Order(2) + 1
                                for i = 1 : NURBS.Order(1) + 1
                                    N0(1, l) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1) * Nz(ez, qz, k, 1);
                                    N1(1, l) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 1) * Nz(ez, qz, k, 1);
                                    N1(2, l) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 2) * Nz(ez, qz, k, 1);
                                    N1(3, l) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1) * Nz(ez, qz, k, 2);
                                    l = l + 1;
                                end
                            end
                        end
                        [R0, R1] = Rationalize(Weights(ElConn), N0, N1);
                        
                        J2 = abs(Jx(ex) * Jy(ey) * Jz(ez));
                        W = Wx(qx) * Wy(qy) * Wz(qz);
                        % Gradient of mapping from parameter space to physical space
                        dxdxi = R1 * CtrlPts(ElConn, :);
                        
                        J1 = abs(det(dxdxi));
                        
                        %Compute derivatives of basis functions w.r.t physical coordinates
                        dRdx = dxdxi^(-1) * R1;
                        
                        B = dRdx; % size(B) = [NSD, NEN]
                        % NEN: number of local basis functions (in fem, n_en = number of element nodes)
                        
                        % compute element stiffness at quadrature point
                        Ke = Ke + B' * D * B * J1 * J2 * W;
                        
                        Pts = R0 * CtrlPts(ElConn, 1 : 3);
                        
                        x = Pts(1); y = Pts(2); z = Pts(3);
                        Fe = Fe + R0' * s(x, y, z) * J1 * J2 * W; % element nodal source vector
                    end
                end
            end          
            FVals(:, e) = Fe;
            KVals(:, e) = Ke(:);
        end
    end
end
end