function [Bstrain, WeJ, countGP]= calcLocalBMatrices2D_mod(Mesh, NURBS, GNum, GDof)
% function [Bstrain, WeJ, countGP] = calcLocalBMatrices2D(Mesh, NURBS)
% -------------------------------------------------------------------
% Calculate two dimensional local B matrices
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
% -------------------------------------------------------------------
% Output:
%       Bstrain: matrix stores local B matrices at all Gauss points
%       (size(Bstrain) = NEl * NGPs
%       WeJ (Nel*NGPs): vector stores J1*J2*W through all Gauss points
%       countGP: total numeber of Gauss points over entire problem
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

assert(iscell(Mesh));

NPatch = numel(Mesh);

NGPs = NURBS{1}.Order + 1; % compatible patches

NEl = 0;
for iPtc = 1 : NPatch
    NEl = NEl + Mesh{iPtc}.NEl;
end
Bstrain = cell(NEl * prod(NGPs), 1);
WeJ = zeros(NEl * prod(NGPs), 1);

countGP = 0;   % count the number of Gauss points

for iPtc = 1 : NPatch
    
    [Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS{iPtc}.Order(1), NURBS{iPtc}.NCtrlPts(1), NURBS{iPtc}.KntVect{1}, 1, NGPs(1), Mesh{iPtc}.NElDir(1));
    [Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS{iPtc}.Order(2), NURBS{iPtc}.NCtrlPts(2), NURBS{iPtc}.KntVect{2}, 1, NGPs(2), Mesh{iPtc}.NElDir(2));
    N0 = zeros(1, Mesh{iPtc}.NEN);
    N1 = zeros(NURBS{iPtc}.Dim, Mesh{iPtc}.NEN);
    
    Weights = reshape(NURBS{iPtc}.Weights, 1, []);
    CtrlPts = reshape(NURBS{iPtc}.CtrlPts3D, 3, [])';
    
    
    for ey = 1 : Mesh{iPtc}.NElDir(2)
        for ex = 1 : Mesh{iPtc}.NElDir(1)
            e = sub2ind(Mesh{iPtc}.NElDir, ex, ey);
            
            connectivity(1 : Mesh{iPtc}.NEN) = GNum{iPtc}(Mesh{iPtc}.El(e, :));
            connectivity(Mesh{iPtc}.NEN + 1 : 2*Mesh{iPtc}.NEN) = GNum{iPtc}(Mesh{iPtc}.El(e, :) + Mesh{iPtc}.NDof / Mesh{iPtc}.Dof);
            
            for qy = 1 : NGPs(2)
                for qx = 1 : NGPs(1)
                    k = 1;
                    for j = 1 : NURBS{iPtc}.Order(2) + 1
                        for i = 1 : NURBS{iPtc}.Order(1) + 1
                            N0(k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1);
                            N1(1, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 1);
                            N1(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 2);
                            k = k + 1;
                        end
                    end
                    [~, R1] = Rationalize(Weights(Mesh{iPtc}.El(e, :)), N0, N1);
                    J2 = Jx(ex) * Jy(ey);
                    W = Wx(qx) * Wy(qy);
                    
                    % Gradient of mapping from parameter space to physical space
                    dxdxi = R1 * CtrlPts(Mesh{iPtc}.El(e, :), 1 : 2);
                    
                    J1 = abs(det(dxdxi));
                    
                    %Compute derivatives of basis functions w.r.t physical coordinates
                    dRdx = dxdxi \ R1;
                    
                    % B matrix
                    B = zeros(3, 2 * Mesh{iPtc}.NEN);
                    
                    % way to sort dispacelement: [u1 u2.....u_Ndof v1,v2,...]
                    B(1, 1 : Mesh{iPtc}.NEN) = dRdx(1, :);
                    B(2, Mesh{iPtc}.NEN + 1 : end) = dRdx(2, :);
                    B(3, 1 : Mesh{iPtc}.NEN) = dRdx(2, :);
                    B(3, Mesh{iPtc}.NEN + 1 : end) = dRdx(1, :);
                    
                    % compute element B and WeJ at quadrature point
                    Btemp = zeros(3, GDof);
                    Btemp(:, connectivity) = B;
                    countGP = countGP + 1;
                    Bstrain{countGP} = sparse(Btemp);
                    WeJ(countGP) = J1 * J2 * W;
                end
            end
        end
    end
end
end
