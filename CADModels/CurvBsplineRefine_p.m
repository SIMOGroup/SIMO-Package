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
% CtrlPts = zeros(2, 5);
% CtrlPts(:, 1) = [0; 0];
% CtrlPts(:, 2) = [0.8; 0];
% CtrlPts(:, 3) = [0.5; 0.5];
% CtrlPts(:, 4) = [1; 1];
% CtrlPts(:, 5) = [1; 0];
CtrlPts = zeros(3, 3);

CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [0.5; 1];
CtrlPts(1 : 2, 3) = [1; 0];

KntVect = [0 0 0.5 1 1]; % knot vector

NCtrlPts = size(CtrlPts, 2);
p = numel(KntVect) - NCtrlPts - 1;
xi = linspace(0, 1, 101);

[KntVect, CtrlPts] = DegreeElevateCurv(p, KntVect, CtrlPts, 0);

NCtrlPts = size(CtrlPts, 2);
p = numel(KntVect) - NCtrlPts - 1;
Idx = FindSpan(NCtrlPts, p, xi, KntVect);
N = BasisFuns(Idx, xi, p, KntVect);
C = zeros(numel(xi), size(CtrlPts, 1));
for ip = 0 : p
    C = C + bsxfun(@times, N(:, ip + 1),...
        CtrlPts(:, Idx - p + ip)');
end

IdxKnt = FindSpan(NCtrlPts, p, unique(KntVect), KntVect);
NKnt = BasisFuns(IdxKnt, unique(KntVect), p, KntVect);
CKnt = zeros(numel(unique(KntVect)), size(CtrlPts, 1));
for ip = 0 : p
    CKnt = CKnt + bsxfun(@times, NKnt(:, ip + 1),...
        CtrlPts(:, IdxKnt - p + ip)');
end

figure
hold on
% Plot curve
plot(C(:, 1), C(:, 2));
% Plot control polygon
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--');
% Plot control points
plot(CtrlPts(1, :), CtrlPts(2, :),...
    'r.','MarkerSize',30);
% Plot 5 points on the curve corresponding to 5 knots
plot(CKnt(:, 1), CKnt(:, 2),'s','MarkerEdgeColor', 'g',...
    'MarkerFaceColor', 'g', 'MarkerSize', 5);

N = zeros(numel(xi), NCtrlPts);
for i = 1 : NCtrlPts
    N(:, i) = OneBasisFun(p, KntVect, i, xi);
end

figure
plot(xi, N);
%axis equal
