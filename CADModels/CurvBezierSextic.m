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

% control points in homogeneous 4D coordinates
CtrlPts = zeros(2, 7);
CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [0; 1];
CtrlPts(1 : 2, 3) = [1; 2];
CtrlPts(1 : 2, 4) = [2.5; -0.5];
CtrlPts(1 : 2, 5) = [4; 2];
CtrlPts(1 : 2, 6) = [5; 2.5];
CtrlPts(1 : 2, 7) = [6; 1];

p = size(CtrlPts, 2) - 1; % degree of the curve

% evalutate the bernstein basis functions
xi = linspace(0, 1, 401);
B = AllBernstein(p, xi);

% contruct bezier curve
C = B * CtrlPts';

% plot basis functions
figure
set(gcf,'color','white')
axis equal
daspect([1 1 1])
set(gcf, 'color', 'white');
plot(xi, B, 'LineWidth', 1.5)

% plot curve, control points and control polygon
figure
hold on
set(gcf, 'color', 'white')
% project to cartesian 3D coordinates
plot(C(:,1), C(:, 2), 'b-', 'LineWidth', 1.5)
plot(CtrlPts(1, :), CtrlPts(2, :),'r.', 'MarkerSize', 15)
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--')
axis equal
axis off

