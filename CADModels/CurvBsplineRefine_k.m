% -------------------------------------------------------------------------
% k-refinement for one dimensional geometry
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

% control points
CtrlPts = zeros(3, 2);

CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [0.5; 1];
CtrlPts(1 : 2, 3) = [1; 0];

KntVect = [0 0 0 1 1 1]; % knot vector

NCtrlPts = size(CtrlPts, 2);
p = numel(KntVect) - NCtrlPts - 1; % order of basis functions

% first, we elevate the order of the coarsest mesh.
[KntVect, CtrlPts] = DegreeElevateCurv(p, KntVect, CtrlPts, 1);

NCtrlPts = size(CtrlPts, 2);
p = numel(KntVect) - NCtrlPts - 1;

% then we insert the knots 
[KntVect, CtrlPts] = RefineKntVectCurv(p, KntVect, ...
    CtrlPts, [1/2]);

NCtrlPts = size(CtrlPts, 2);
p = numel(KntVect) - NCtrlPts - 1;

xi = linspace(0, 1, 301); % parametric points
% Evaluate the basis functions
N = zeros(numel(xi), NCtrlPts);
for i = 1 : NCtrlPts
    N(:, i) = OneBasisFun(p, KntVect, i, xi);
end

figure
%hold on
%set(gcf,'color','white')
% set(gca,'XTick', 0:1/3:1)
% set(gca,'YTick', 0:1:1)
% daspect([1 2 1])
plot(xi, N);
%axis tight


