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

% radius
radius = 1;

CtrlPts = zeros(4, 6);
CtrlPts(1 : 3, 1) = [0; -1; 0];
CtrlPts(1 : 3, 2) = [-4; -1; 0];
CtrlPts(1 : 3, 3) = [-2; 3; 0];
CtrlPts(1 : 3, 4) = [2; 3; 0];
CtrlPts(1 : 3, 5) = [4; -1; 0];
CtrlPts(1 : 3, 6) = [0; -1; 0];

CtrlPts(1 : 3, :) = CtrlPts(1 : 3, :) * radius;

% weights
CtrlPts(4, :) = 1;
fac = 5;
CtrlPts(:, 1 : 5 : 6) = CtrlPts(:, 1 : 5 : 6) * fac;

knots{1} = [0 0 0 0 0 0 1 1 1 1 1 1];

curv = CreateNURBS(knots, CtrlPts);

figure
hold on
axis off
set(gcf, 'color', 'white');
daspect([1 1 1]);
axis equal
PlotGeo(curv)
PlotKnts(curv)
PlotCtrlPts(curv)
PlotCtrlNet(curv)

