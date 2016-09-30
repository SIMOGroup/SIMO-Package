% -------------------------------------------------------------------------
% An example of a quintic Bezier curve as circle
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

% radius of the circle
radius = 4;

% control points in homogeneous 4D coordinates
CtrlPts = zeros(4, 6);
CtrlPts(1 : 3, 1) = [0; -1; 0];
CtrlPts(1 : 3, 2) = [-4; -1; 0];
CtrlPts(1 : 3, 3) = [-2; 3; 0];
CtrlPts(1 : 3, 4) = [2; 3; 0];
CtrlPts(1 : 3, 5) = [4; -1; 0];
CtrlPts(1 : 3, 6) = [0; -1; 0];

CtrlPts(1 : 3, :) = CtrlPts(1 : 3, :) * radius;
% assign the weights
CtrlPts(4, :) = 1;
fac = 5;
CtrlPts(:, 1 : 5 : 6) = CtrlPts(:, 1 : 5 : 6) * fac;

p = size(CtrlPts, 2) - 1;%degree of the curve

% evalutate the bernstein basis functions
xi = linspace(0, 1, 101);
B = AllBernstein(p, xi);

% contruct the rational bezier curve in 4D space
Cw = CtrlPts * B';

% project it to 3D space
C = bsxfun(@rdivide, Cw, Cw(4, :));

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
plot(C(1, :), C(2, :), 'b-', 'LineWidth', 1.5)
weights = CtrlPts(4, :);
plot(CtrlPts(1, :)./weights, CtrlPts(2, :)./weights,...
    'r.', 'MarkerSize', 15)
plot(CtrlPts(1, :)./weights, CtrlPts(2, :)./weights, 'k--')
axis equal

