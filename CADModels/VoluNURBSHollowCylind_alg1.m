% -------------------------------------------------------------------------
% representing model of the hollow cylinder
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

r = 2; % inner radius of hollow cylinder
t = 0.5; % thickness of hollow cylinder
h = 6; % height of hollow cylinder

% control points
CtrlPts = zeros(4, 9, 2, 2);

CtrlPts(1 : 3, 1, 1, 1) = [ r;  0; 0];
CtrlPts(1 : 3, 2, 1, 1) = [ r;  r; 0];
CtrlPts(1 : 3, 3, 1, 1) = [ 0;  r; 0];
CtrlPts(1 : 3, 4, 1, 1) = [-r;  r; 0];
CtrlPts(1 : 3, 5, 1, 1) = [-r;  0; 0];
CtrlPts(1 : 3, 6, 1, 1) = [-r; -r; 0];
CtrlPts(1 : 3, 7, 1, 1) = [ 0; -r; 0];
CtrlPts(1 : 3, 8, 1, 1) = [ r; -r; 0];
CtrlPts(1 : 3, 9, 1, 1) = [ r;  0; 0];

CtrlPts(1 : 3, 1, 2, 1) = [ r + t;      0; 0];
CtrlPts(1 : 3, 2, 2, 1) = [ r + t;  r + t; 0];
CtrlPts(1 : 3, 3, 2, 1) = [     0;  r + t; 0];
CtrlPts(1 : 3, 4, 2, 1) = [-(r + t);    r + t; 0];
CtrlPts(1 : 3, 5, 2, 1) = [-(r + t);        0; 0];
CtrlPts(1 : 3, 6, 2, 1) = [-(r + t); -(r + t); 0];
CtrlPts(1 : 3, 7, 2, 1) = [     0; -(r + t); 0];
CtrlPts(1 : 3, 8, 2, 1) = [ r + t; -(r + t); 0];
CtrlPts(1 : 3, 9, 2, 1) = [ r + t;        0; 0];

CtrlPts(1 : 3, 1, 1, 2) = [ r;  0; h];
CtrlPts(1 : 3, 2, 1, 2) = [ r;  r; h];
CtrlPts(1 : 3, 3, 1, 2) = [ 0;  r; h];
CtrlPts(1 : 3, 4, 1, 2) = [-r;  r; h];
CtrlPts(1 : 3, 5, 1, 2) = [-r;  0; h];
CtrlPts(1 : 3, 6, 1, 2) = [-r; -r; h];
CtrlPts(1 : 3, 7, 1, 2) = [ 0; -r; h];
CtrlPts(1 : 3, 8, 1, 2) = [ r; -r; h];
CtrlPts(1 : 3, 9, 1, 2) = [ r;  0; h];

CtrlPts(1 : 3, 1, 2, 2) = [ r + t;      0; h];
CtrlPts(1 : 3, 2, 2, 2) = [ r + t;  r + t; h];
CtrlPts(1 : 3, 3, 2, 2) = [     0;  r + t; h];
CtrlPts(1 : 3, 4, 2, 2) = [-(r + t);    r + t; h];
CtrlPts(1 : 3, 5, 2, 2) = [-(r + t);        0; h];
CtrlPts(1 : 3, 6, 2, 2) = [-(r + t); -(r + t); h];
CtrlPts(1 : 3, 7, 2, 2) = [     0; -(r + t); h];
CtrlPts(1 : 3, 8, 2, 2) = [ r + t; -(r + t); h];
CtrlPts(1 : 3, 9, 2, 2) = [ r + t;        0; h];
% weights
CtrlPts(4, :, :, :) = 1;
fac = 1 / sqrt(2);
CtrlPts(:, 2 : 2 : 8, :, :) = CtrlPts(:, 2 : 2 : 8, :, :) * fac;

% knot vector
KntVect{1} = [0 0 0 1 1 2 2 3 3 4 4 4];
KntVect{1} = KntVect{1}./max(KntVect{1});
KntVect{2} = [0 0 1 1];
KntVect{3} = [0 0 1 1];

Volu = CreateNURBS(KntVect, CtrlPts);

figure
hold on
axis off
view(3)
daspect([1, 1, 1])
rotate3d on
set(gcf,'color','white')

PlotGeo(Volu)
PlotCtrlPts(Volu)
PlotKnts(Volu)
PlotCtrlNet(Volu)
