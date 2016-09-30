function KVals = calcLocalStiffnessMatrices3D(Mesh, NURBS, E, nu)
% function KVals = calcLocalStiffnessMatrices3D(Mesh, NURBS, E, nu)
% -------------------------------------------------------------------
% Calculate three dimensional stiffness matrices
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       E: Young's modulus
%       nu: Poisson's ratio
% -------------------------------------------------------------------
% Output:
%       KVals: matrix stores local stiffness matrices
%           (size(KVals) = [(NEN * 3) ^ 2, NEL])
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

D = getElastMat(E, nu, '3D');
NGPs = NURBS.Order + 1;
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));
[Jz, Wz, ~, Nz] = calcDersBasisFunsAtGPs(NURBS.Order(3), NURBS.NCtrlPts(3), NURBS.KntVect{3}, 1, NGPs(3), Mesh.NElDir(3));

N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);

KVals = zeros((Mesh.NEN * Mesh.Dof) ^ 2, Mesh.NEl);

CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';
Weights = reshape(NURBS.Weights, 1, []);

for ez = 1 : Mesh.NElDir(3)
    for ey = 1 : Mesh.NElDir(2)
        for ex = 1 : Mesh.NElDir(1)
            e = sub2ind(Mesh.NElDir, ex, ey, ez);
            Ke = zeros(Mesh.NEN * Mesh.Dof);
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
                        [~, R1] = Rationalize(Weights(Mesh.El(e, :)), N0, N1);
                        
                        J2 = Jx(ex) * Jy(ey) * Jz(ez);
                        W = Wx(qx) * Wy(qy) * Wz(qz);
                        % Gradient of mapping from parameter space to physical space
                        dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
                        
                        J1 = abs(det(dxdxi));
                        
                        %Compute derivatives of basis functions w.r.t physical coordinates
                        dRdx = dxdxi^(-1) * R1;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % COMPUTE B MATRIX
                        % %       _                                                            _
                        % %     | N_{1, x} N_{2, x} ... 0        0        ... 0        0        ...|
                        % %     | 0        0        ... N_{1, y} N_{2, y} ... 0        0        ...|
                        % % B = | 0        0        ... 0        0        ... N_{1, z} N_{2, z} ...|
                        % %     | N_{1, y} N_{2, y} ... N_{1, x} N_{2, x} ... 0        0        ...|
                        % %     | 0        0        ... N_{1, z} N_{2, z} ... N_{1, y} N_{2, y} ...|
                        % %     | N_{1, z} N_{2, z} ... 0        0        ... N_{1, x} N_{2, x} ...|
                        % %       -                                                            -
                        % % \sigma = [\sigma_{xx} \sigma_{yy} \sigma_{zz} \sigma_{xy} \sigma_{yz} \sigma_{xz}]
                        B = zeros(6, 3 * Mesh.NEN);
                        B(1, 1 : Mesh.NEN) = dRdx(1, :);
                        B(2,  Mesh.NEN + 1 : 2 * Mesh.NEN) = dRdx(2, :);
                        B(3, 2 * Mesh.NEN + 1 : end) = dRdx(3, :);
                        
                        B(4, 1 : Mesh.NEN) = dRdx(2, :);
                        B(4, Mesh.NEN + 1 : 2 * Mesh.NEN) = dRdx(1, :);
                        
                        B(5, 2 * Mesh.NEN + 1 : end)  = dRdx(2, :);
                        B(5, Mesh.NEN + 1 : 2 * Mesh.NEN) = dRdx(3, :);
                        
                        B(6, 1 : Mesh.NEN) = dRdx(3, :);
                        B(6, 2 * Mesh.NEN + 1 : end) = dRdx(1, :);
                        
                        % compute element stiffness at quadrature point
                        Ke = Ke + B' * D * B * J1 * J2 * W;
                    end
                end
            end
            KVals(:, e) = Ke(:);
        end
    end
end
end