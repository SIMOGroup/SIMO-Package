% -------------------------------------------------------------------------
% Plot the B-spline curve
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

close all
clear
clc

% control points
% CtrlPts = zeros(2, 7);
% CtrlPts(1 : 2, 1) = [0; 0];
% CtrlPts(1 : 2, 2) = [0; 1];
% CtrlPts(1 : 2, 3) = [1; 1];
% CtrlPts(1 : 2, 4) = [2.5; -0.5];
% CtrlPts(1 : 2, 5) = [4; 2];
% CtrlPts(1 : 2, 6) = [5; 2.5];
% CtrlPts(1 : 2, 7) = [6; 1];

%CtrlPts = zeros(2, 4);
CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [1/4; 1/2];
CtrlPts(1 : 2, 3) = [3/4; 3/4];
CtrlPts(1 : 2, 4) = [1; 0];

% knot vector
%KntVect = [0 0 0 0 1 2 3 4 4 4 4];
KntVect = [0 0 0 1/2 1 1 1];
KntVect = KntVect ./ max(KntVect);

NCtrlPts = size(CtrlPts, 2); % number of control points
p = numel(KntVect) - NCtrlPts - 1; % order of basis functions
ParaPts = linspace(0, 1, 101); % paramatric points

Idx = FindSpan(NCtrlPts, p, ParaPts, KntVect);
N0 = BasisFuns(Idx, ParaPts, p, KntVect);

C = zeros(2, numel(ParaPts));
for i = 1 : p + 1
    C = C + bsxfun(@times, N0(:, i)', CtrlPts(:, Idx - p + i - 1));
end

figure
hold on
daspect([1 1 1])
axis equal
axis off
% Plot curve
plot(C(1, :), C(2, :));
% Plot control polygon
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--');
% Plot control points
plot(CtrlPts(1, :), CtrlPts(2, :),...
    'r.','MarkerSize',15);
