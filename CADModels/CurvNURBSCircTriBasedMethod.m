%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% representing a full circle using triangle based method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

radius = 1;

% Control points in homogeneous space

CtrlPts = zeros(4, 7);
CtrlPts(1 : 2, 1, 1) = [cos(pi / 6); sin(pi / 6)];
CtrlPts(1 : 2, 2, 1) = [0; 1 / cos(pi / 3)];
CtrlPts(1 : 2, 3, 1) = [-cos(pi / 6); sin(pi / 6)];
CtrlPts(1 : 2, 4, 1) = [-tan(pi / 3); -1];
CtrlPts(1 : 2, 5, 1) = [0; -1];
CtrlPts(1 : 2, 6, 1) = [tan(pi / 3); -1];
CtrlPts(1 : 2, 7, 1) = [cos(pi / 6); sin(pi / 6)];

CtrlPts(1 : 3, :, :) = CtrlPts(1 : 3, :, :) * radius;
% assign the weights
CtrlPts(4, :, :) = 1;
W = 1 / 2;
CtrlPts(:, 2 : 2 : end, :) = CtrlPts(:, 2 : 2 : end, :) * W;

% knot vector
KntVect{1} = [0 0 0 1/3 1/3 2/3 2/3 1 1 1];

Curv = CreateNURBS(KntVect, CtrlPts);

figure
hold on
axis off
set(gcf, 'color', 'white');
daspect([1 1 1]);
axis equal
PlotGeo(Curv)
PlotKnts(Curv)
PlotCtrlPts(Curv)
PlotCtrlNet(Curv)
