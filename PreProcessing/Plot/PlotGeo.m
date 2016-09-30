function PlotGeo(NURBS)
% PlotGeo(NURBS)
% Plot geometry
% ------------------------------------------------------------

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

if NURBS.Dim == 1
    xi = linspace(0, 1, 41); % parametric points
    Cw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, {xi});
    w = Cw(4, :);
    C = bsxfun(@rdivide, Cw, w);
    plot3 (C(1, :), C(2, :), C(3, :), 'Linewidth', 1, 'color', 'k');
elseif NURBS.Dim == 2
    xi = linspace(0, 1, 101); % parametric points
    eta = linspace(0, 1, 101); % parametric points
    % elvaluate the parametric points to plot the surface
    Sw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, {xi, eta});
    % ----------------------------------------------------------
    
    % Plot the NURBS surface
    bcol=[173 234 234] ./ 250;
    bcol = repmat(bcol, 3, 1);
    colormap(bcol);
%     colormap('white');
%     light = 'on';
    [~, m, n] = size(Sw);
    w = Sw(4, :, :);
    % parametric points projected into Cartesian 3D space.
    S = bsxfun(@rdivide, Sw, w);
    x = reshape(S(1, :, :), m, n);
    y = reshape(S(2, :, :), m, n);
    z = reshape(S(3, :, :), m, n);
    surfl(x, y, z);
    shading interp
elseif NURBS.Dim == 3
    for dir = 1 : 3
        for side = [0 1]
            Surf = NURBSExtract(NURBS, dir, side);
            PlotGeo(Surf);
        end
    end    
end
end