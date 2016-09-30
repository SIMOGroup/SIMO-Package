% -------------------------------------------------------------------------
% representing a full circle using triangle based method
% -------------------------------------------------------------------------

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

clear
close all
clc

radius = 1;

% Control points in homogeneous space

CtrlPts = zeros(4, 7);
CtrlPts(1 : 2, 1, 1) = [cos(pi / 6); sin(pi / 6)];
CtrlPts(1 : 2, 2, 1) = [0; 1 / cos(pi / 3)];
CtrlPts(1 : 2, 3, 1) = [-cos(pi / 6); sin(pi / 6)];
CtrlPts(1 : 2, 4, 1) = [-tan(pi / 3); -1];
CtrlPts(1 : 2, 5, 1) = [0; -1];
CtrlPts(1 : 2, 6, 1) = [tan(pi / 3); -1];
CtrlPts(1 : 2, 7, 1) = [cos(pi / 6); sin(pi / 6)];

CtrlPts(1 : 3, :, :) = CtrlPts(1 : 3, :, :) * radius;
% assign the weights
CtrlPts(4, :, :) = 1;
W = 1 / 2;
CtrlPts(:, 2 : 2 : end, :) = CtrlPts(:, 2 : 2 : end, :) * W;

% knot vector
KntVect{1} = [0 0 0 1/3 1/3 2/3 2/3 1 1 1];

Curv = CreateNURBS(KntVect, CtrlPts);

figure
hold on
axis off
set(gcf, 'color', 'white');
daspect([1 1 1]);
axis equal
PlotGeo(Curv)
PlotKnts(Curv)
PlotCtrlPts(Curv)
PlotCtrlNet(Curv)
