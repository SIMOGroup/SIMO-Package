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

clc
clear
close all

% control points

radius = 12.7 / 2;

CtrlPts = zeros(4, 3, 3);
CtrlPts(1 : 2, 1, 1) = [-cos(pi / 4); sin(pi / 4)];
CtrlPts(1 : 2, 2, 1) = [-1 / cos(pi / 4); 0];
CtrlPts(1 : 2, 3, 1) = [-cos(pi / 4); -sin(pi / 4)];

CtrlPts(1 : 2, 1, 2) = [0; 1 / cos(pi / 4)];
CtrlPts(1 : 2, 2, 2) = [0; 0];
CtrlPts(1 : 2, 3, 2) = [0; -1 / cos(pi/4)];

CtrlPts(1 : 2, 1, 3) = [cos(pi / 4); sin(pi / 4)];
CtrlPts(1 : 2, 2, 3) = [1 / cos(pi / 4); 0];
CtrlPts(1 : 2, 3, 3) = [cos(pi/4); -sin(pi/4)];

CtrlPts(1 : 2, :, :) = CtrlPts(1 : 2, :, :) * radius;
% weights
CtrlPts(4, :, :) = 1;
W = 1 / sqrt(2);
CtrlPts(:, 2, 1 : 2 : 3) = CtrlPts(:, 2, 1 : 2 : 3) * W;
CtrlPts(:, 1 : 2 : 3, 2) = CtrlPts(:, 1 : 2 : 3, 2) * W;

% knot vectors
KntVect{1} = [0 0 0 1 1 1];
KntVect{2} = [0 0 0 1 1 1];

Surf = CreateNURBS(KntVect, CtrlPts);

% generatrix
CtrlPts = zeros(4, 6);

CtrlPts(1 : 3, 1) = [0; 6.35; 0];
CtrlPts(1 : 3, 2) = [0; 3.18; 7.16];
CtrlPts(1 : 3, 3) = [0; 3.18; 15];
CtrlPts(1 : 3, 4) = [0; 3.18; 45];
CtrlPts(1 : 3, 5) = [0; 3.18; 52.84];
CtrlPts(1 : 3, 6) = [0; 6.35; 60];

% weights
CtrlPts(4, :) = 1;
Ws = cos(deg2rad(23.9));
CtrlPts(:, 2 : 3 : 5) = CtrlPts(:, 2 : 3 : 5) * Ws;

KntVect1{1} = [0 0 0 1 2 3 4 4 4];
KntVect1{1} = KntVect1{1} ./ max(KntVect1{1});

Curv = CreateNURBS(KntVect1, CtrlPts);

volu = NURBSSweep(Surf, Curv);

figure
hold on
grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
PlotGeo(volu)
PlotKnts(volu)
PlotCtrlPts(volu)
PlotCtrlNet(volu)
