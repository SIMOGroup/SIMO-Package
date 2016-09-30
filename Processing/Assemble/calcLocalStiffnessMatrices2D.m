function KVals = calcLocalStiffnessMatrices2D(Mesh, NURBS, E, nu, lab)
% function KVals = calcLocalStiffnessMatrices2D(Mesh, NURBS, E, nu, lab)
% -------------------------------------------------------------------
% Calculate two dimensional local stiffness matrices
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       E: Young's modulus
%       nu: Poisson's ratio
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

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));
N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);

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
                        k = k + 1;
                    end
                end
                [~, R1] = Rationalize(Weights(Mesh.El(e, :)), N0, N1);
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % Gradient of mapping from parameter space to physical space
                dxdxi = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
                
                J1 = abs(det(dxdxi));
                
                %Compute derivatives of basis functions w.r.t physical coordinates
                dRdx = dxdxi^(-1) * R1;
                
                % B matrix
                %        _                                      _
                %        |  N_1,x  N_2,x  ...      0      0  ... |
                %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
                %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
                %        -                                      -
                B = zeros(3, 2 * Mesh.NEN);
                B(1, 1 : Mesh.NEN) = dRdx(1, :);
                B(2, Mesh.NEN + 1 : end) = dRdx(2, :);
                B(3, 1 : Mesh.NEN) = dRdx(2, :);
                B(3, Mesh.NEN + 1 : end) = dRdx(1, :);
                
                % compute element stiffness at quadrature point
                Ke = Ke + B' * D * B * J1 * J2 * W;
            end
        end
        KVals(:, e) = Ke(:);
    end
end
end