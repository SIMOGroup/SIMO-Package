function PlotFlux2D(NURBS, Mesh, d, ka)
% function PlotFlux2D(NURBS, Mesh, d, ka)
% Plot heat flux at gauss points

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

NGPs = NURBS.Order + 1;

[~, ~, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[~, ~, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));

N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);

D  = ka * eye(2);

Weights = reshape(NURBS.Weights, 1, []);
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';

pts = zeros(prod(NGPs), 2);
q = size(pts);

figure
hold on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
axis off
PlotGeo(NURBS);
PlotKnts(NURBS);
for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        e = sub2ind([Mesh.NElDir(1), Mesh.NElDir(2)], ex, ey);
        ElConn = Mesh.El(e, :); %  element connectivity        
        de = d(ElConn);     
        % compute flux vector
        ind = 1;
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
                [R0, R1] = Rationalize(Weights(ElConn), N0, N1);
                
                % Gradient of mapping from parameter space to physical space
                dxdxi = R1 * CtrlPts(ElConn, 1 : 2);
                
                % Compute derivatives of basis functions w.r.t physical coordinates
                dRdx = dxdxi^(-1) * R1;
                
                B = dRdx;
                
                % gauss points in physical coordinates
                pts(ind, :) = R0 * CtrlPts(ElConn, 1 : 2);
                
                q(ind, :) = -D * B * de; % compute the flux
                ind = ind + 1;
            end
        end
        quiver(pts(:, 1), pts(:, 2), q(:, 1), q(:, 2), 'k');
%         plot(pts(:, 1), pts(:, 2), 'k*');
        title('Heat Flux');
        xlabel('X');
        ylabel('Y');
    end
end

end