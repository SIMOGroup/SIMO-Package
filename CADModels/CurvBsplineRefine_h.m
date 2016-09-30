
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
% CtrlPts = zeros(4, 5);
% CtrlPts(1 : 2, 1) = [0; 0];
% CtrlPts(1 : 2, 2) = [0.8; 0];
% CtrlPts(1 : 2, 3) = [0.5; 0.5];
% CtrlPts(1 : 2, 4) = [1; 1];
% CtrlPts(1 : 2, 5) = [1; 0];
% CtrlPts(4, :) = 1;
% KntVect{1} = [0 0 0 0 0 1 1 1 1 1]; % knot vector


CtrlPts = zeros(2, 2);

CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [1; 0];

CtrlPts(4, :) = 1;

KntVect{1} = [0 0  1 1]; % knot vector

NCtrlPts = size(CtrlPts, 2);
p = numel(KntVect{1}) - NCtrlPts - 1;
nel=2;
insertKnot=[1/nel:1/nel:1-1/nel];

[refineKnts, finCtrlPts] = RefineKntVectCurv(p, KntVect{1}, CtrlPts,insertKnot);

curv =  CreateNURBS({refineKnts}, finCtrlPts);
figure
hold on
axis off
set(gcf,'color','white')

PlotGeo(curv)
PlotCtrlPts(curv)
PlotKnts(curv)
PlotCtrlNet(curv)


% plot basis functions

NCtrlPts = size(finCtrlPts, 2);
p = numel(refineKnts) - NCtrlPts - 1;
xi = linspace(0, 1, 101);

N0 = zeros(numel(xi), NCtrlPts);
for i = 1 : NCtrlPts
    N0(:, i) = OneBasisFun(p, refineKnts, i, xi);
end

figure
plot(xi, N0);
%
