% -------------------------------------------------------------------------
% Plot a cubic B-spline curve
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

% controlpoints
CtrlPts = zeros(2, 7);
CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [0; 1];
CtrlPts(1 : 2, 3) = [1; 2];
CtrlPts(1 : 2, 4) = [2.5; -0.5];
CtrlPts(1 : 2, 5) = [4; 2];
CtrlPts(1 : 2, 6) = [5; 2.5];
CtrlPts(1 : 2, 7) = [6; 1];

% knot vector
KntVect = [0 0 0 0 1 2 3 4 4 4 4];
KntVect = KntVect ./ max(KntVect);

xi = linspace(0, 1, 101); % parametric points

% elvaluate the parametric points to plot the curve
C = BsplineEval({KntVect}, CtrlPts, {xi});

% elvaluate the knots
uqKntVect = KntVect([true; diff(KntVect(:)) > 0]); 
evalKnts = BsplineEval({KntVect}, CtrlPts, {uqKntVect});

figure
hold on
set(gcf, 'color', 'white');
axis off
axis equal
% Plot curve
plot(C(1, :), C(2, :), 'linewidth', 1.5);
% Plot control polygon
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--');
% Plot control points
plot(CtrlPts(1, :), CtrlPts(2, :),...
    'r.','MarkerSize',15);
% Plot 5 points on the curve corresponding to 5 knots
plot(evalKnts(1, :), evalKnts(2, :),'s','MarkerEdgeColor', 'g',...
    'MarkerFaceColor', 'g', 'MarkerSize', 5);
axis equal
