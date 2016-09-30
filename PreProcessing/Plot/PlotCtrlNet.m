function PlotCtrlNet(NURBS)
% function PlotCtrlNet(NURBS)
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
    CtrlPts3D = NURBS.CtrlPts3D(1 : 3, :);
    plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
elseif NURBS.Dim == 2
    % Plot control points along y-direction
    for i = 1 : NURBS.NCtrlPts(1)
        CtrlPts3D = reshape(NURBS.CtrlPts3D(:, i, :), 3, []);
        plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--');
    end
    
    % Plot control points along x-direction
    for j = 1 : NURBS.NCtrlPts(2)
        CtrlPts3D = reshape(NURBS.CtrlPts3D(:, :, j), 3, []);
        plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
    end
elseif NURBS.Dim == 3
    for i = 1 : NURBS.NCtrlPts(1)
        for j = 1 : NURBS.NCtrlPts(2)
            CtrlPts3D = reshape(NURBS.CtrlPts3D(:, i, j, :), 3, []);
            plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
        end
        for k = 1 : NURBS.NCtrlPts(3)
            CtrlPts3D = reshape(NURBS.CtrlPts3D(:, i, :, k), 3, []);
            plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
        end
    end
    for j = 1 : NURBS.NCtrlPts(2)
        for k = 1 : NURBS.NCtrlPts(3)
            CtrlPts3D = reshape(NURBS.CtrlPts3D(:, :, j, k), 3, []);
            plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
        end
    end    
end
end