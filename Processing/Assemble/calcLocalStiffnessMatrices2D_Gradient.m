function KVals = calcLocalStiffnessMatrices2D_Gradient(Mesh, NURBS, E, nu, l, lab)
% function KVals = calcLocalStiffnessMatrices2D_Gradient(Mesh, NURBS, E, nu, l, lab)
% -------------------------------------------------------------------
% Calculate two dimensional local stiffness matrices for gradient
% elasticity problems
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       E: Young's modulus
%       nu: Poisson's ratio
%       l: internal length
%       lab: a string to determine stress state
%             ex. 'PlaneStress', 'PlaneStrain'
% -------------------------------------------------------------------
% Output:
%       KVals: matrix stores local stiffness matrices
%           (size(KVals) = [(NEN * 2) ^ 2, NEL])
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

D = getElastMat(E, nu, lab);

NGPs = NURBS.Order + 1;

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 2, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 2, NGPs(2), Mesh.NElDir(2));
N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);
N2 = zeros(3, Mesh.NEN);

KVals = zeros((Mesh.NEN * Mesh.Dof) ^ 2, Mesh.NEl);

Weights = reshape(NURBS.Weights, 1, []);
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';

for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        e = sub2ind(Mesh.NElDir, ex, ey);
        Ke = zeros(Mesh.NEN * 2);
        for qy = 1 : NGPs(2)
            for qx = 1 : NGPs(1)
                k = 1;
                for j = 1 : NURBS.Order(2) + 1
                    for i = 1 : NURBS.Order(1) + 1
                        N0(k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1);
                        N1(1, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 1);
                        N1(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 2);
                        N2(1, k) = Nx(ex, qx, i, 3) * Ny(ey, qy, j, 1); % d2N/dxi2
                        N2(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 3); % d2N/deta2
                        N2(3, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 2); % d2N/dxi deta
                        k = k + 1;
                    end
                end
                [R0, R1, R2] = Rationalize(Weights(Mesh.El(e, :)), N0, N1, N2);
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % Gradient of mapping from parameter space to physical space
                jacob  = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
                jacob2 = R2 * CtrlPts(Mesh.El(e, :), 1 : 2);
                
                J1 = abs(det(jacob));
                
                dxdxi = jacob(1, 1); 
                dydxi = jacob(1, 2);
                dxdeta = jacob(2, 1); 
                dydeta = jacob(2, 2);
                
                j33 = [dxdxi^2, dydxi^2, 2*dxdxi*dydxi;
                    dxdeta^2, dydeta^2, 2*dxdeta*dydeta;
                    dxdxi*dxdeta, dydxi*dydeta, dxdxi*dydeta + dxdeta*dydxi];
                
                % Jacobian inverse and spatial 1st and 2nd derivatives
                % Compute derivatives of basis functions w.r.t
                % physical coordinates
                dRdx = jacob \ R1;
                dR2dx = j33 \ (R2 - jacob2*dRdx);
                
                % B matrix
                %        _                                      _
                %        |  N_1,x  N_2,x  ...      0      0  ... |
                %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
                %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
                %        -                                      -
                B = zeros(3, 2*Mesh.NEN);
                B(1, 1 : Mesh.NEN) = dRdx(1, :);
                B(2, Mesh.NEN + 1 : end) = dRdx(2, :);
                B(3, 1 : Mesh.NEN) = dRdx(2, :);
                B(3, Mesh.NEN + 1 : end) = dRdx(1, :);
                
                dBdx = zeros(3, 2 * Mesh.NEN, 2);
                dBdx(1, 1 : Mesh.NEN, 1) = dR2dx(1, :);
                dBdx(2, Mesh.NEN + 1 : end, 1) = dR2dx(3, :);
                dBdx(3, 1 : Mesh.NEN, 1) = dR2dx(3, :);
                dBdx(3, Mesh.NEN + 1 : end, 1) = dR2dx(1, :);
                
                dBdx(1, 1 : Mesh.NEN, 2) = dR2dx(3, :);
                dBdx(2, Mesh.NEN + 1 : end, 2) = dR2dx(2, :);
                dBdx(3, 1 : Mesh.NEN, 2) = dR2dx(2, :);
                dBdx(3, Mesh.NEN + 1 : end, 2) = dR2dx(3, :);
                
                ke2 = zeros(Mesh.NEN * 2, Mesh.NEN * 2, 2);
                for i = 1 : 2
                    ke2(:, :, i) = dBdx(:, :, i)'*l^2*D*dBdx(:, :, i)*J1*J2*W;
                end
                % compute element stiffness at quadrature point
                Ke = Ke + B' * D * B * J1 * J2 * W + sum(ke2, 3);
            end
        end
        KVals(:, e) = Ke(:);
    end
end
end