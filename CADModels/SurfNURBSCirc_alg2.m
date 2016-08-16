
clear
close all
clc

% control points

radius = 3;

CtrlPts = zeros(4, 3, 3);
CtrlPts(1 : 2, 1, 1) = [0; 1];
CtrlPts(1 : 2, 2, 1) = [-1; 1];
CtrlPts(1 : 2, 3, 1) = [-1; 0];

CtrlPts(1 : 2, 1, 2) = [1; 1];
CtrlPts(1 : 2, 2, 2) = [0; 0];
CtrlPts(1 : 2, 3, 2) = [-1; -1];

CtrlPts(1 : 2, 1, 3) = [1; 0];
CtrlPts(1 : 2, 2, 3) = [1; -1];
CtrlPts(1 : 2, 3, 3) = [0; -1];

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
%circ_surf = k_refine(Surf, [10 10], [2 2], [1 1]);

figure
hold on
grid on
set(gcf,'color','white')
daspect([1 1 1])
axis([-radius - 0.2 radius + 0.2 -radius - 0.2 radius + 0.2]);
axis equal
% axis off
PlotGeo(Surf)
PlotKnts(Surf)
PlotCtrlPts(Surf)
PlotCtrlNet(Surf)
