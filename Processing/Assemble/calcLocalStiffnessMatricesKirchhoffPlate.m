function [KVals, FVals] = calcLocalStiffnessMatricesKirchhoffPlate(Mesh, NURBS, t, E, nu, q)
% [KVals, FVals] = calcLocalStiffnessMatricesKirchhoffPlate(Mesh, NURBS, t, E, nu, q)
% --------------------------------------------------------------------
% Evaluate local stiffnes matrices and distributed forces for Kirchhoff
% Plate
% --------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       t: thickness of the plate
%       E: Young's modulus
%       nu: Poisson's ratio
%       q: distributed force
% --------------------------------------------------------------------
% Output:
%       KVals: matrix stores local stiffness matrices
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

% constitutive matrix

D  = E * t ^ 3 / (12 * (1 - nu ^ 2));
C  = D * [1 nu 0; nu 1 0; 0 0 0.5 * (1 - nu)];

NGPs = NURBS.Order + 1;

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 2, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 2, NGPs(2), Mesh.NElDir(2));

N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);
N2 = zeros(3, Mesh.NEN);

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
                        N2(1, k) = Nx(ex, qx, i, 3) * Ny(ey, qy, j, 1); % d2N/dxi2
                        N2(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 3); % d2N/deta2
                        N2(3, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 2); % d2N/dxi deta
                        k = k + 1;
                    end
                end
                [R0, R1, R2] = Rationalize(Weights(Mesh.El(e, :)), N0, N1, N2);
                
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % Gradient of mapping from parameter space to physical
                % space
                jacob  = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
                jacob2 = R2 * CtrlPts(Mesh.El(e, :), 1 : 2);
                
                J1    = det(jacob);
                
                dxdxi = jacob(1,1); dydxi = jacob(1,2);
                dxdet = jacob(2,1); dydet = jacob(2,2);
                
                j33   = [dxdxi^2     dydxi^2     2*dxdxi*dydxi;
                    dxdet^2     dydet^2     2*dxdet*dydet;
                    dxdxi*dxdet dydxi*dydet dxdxi*dydet+dxdet*dydxi];
                
                % Jacobian inverse and spatial 1st and 2nd derivatives
                % Compute derivatives of basis functions w.r.t
                % physical coordinates
                dRdx       = jacob \ R1;
                dR2dx      = j33 \ (R2 - jacob2 * dRdx);
                
                % B matrix
                
                B          = dR2dx;
                B(3,:)     = B(3,:)*2;
                
                % compute element stiffness at quadrature point
                Ke = Ke + B' * C * B * J1 * J2 * W;
                
%                 Pts = R0 * CtrlPts(Mesh.El(e, :), :);
%                 
%                 x = Pts(1); y = Pts(2);
                
                % element source vector
                Fe = Fe + q * R0' * J1 * J2 * W;
            end
        end
        FVals(:, e) = Fe;
        KVals(:, e) = Ke(:);
    end
end
end