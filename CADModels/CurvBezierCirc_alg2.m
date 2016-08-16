clear
close all
clc

% radius
radius = 1;

CtrlPts = zeros(4, 6);
CtrlPts(1 : 3, 1) = [0; -1; 0];
CtrlPts(1 : 3, 2) = [-4; -1; 0];
CtrlPts(1 : 3, 3) = [-2; 3; 0];
CtrlPts(1 : 3, 4) = [2; 3; 0];
CtrlPts(1 : 3, 5) = [4; -1; 0];
CtrlPts(1 : 3, 6) = [0; -1; 0];

CtrlPts(1 : 3, :) = CtrlPts(1 : 3, :) * radius;

% weights
CtrlPts(4, :) = 1;
fac = 5;
CtrlPts(:, 1 : 5 : 6) = CtrlPts(:, 1 : 5 : 6) * fac;

knots{1} = [0 0 0 0 0 0 1 1 1 1 1 1];

curv = CreateNURBS(knots, CtrlPts);

figure
hold on
axis off
set(gcf, 'color', 'white');
daspect([1 1 1]);
axis equal
PlotGeo(curv)
PlotKnts(curv)
PlotCtrlPts(curv)
PlotCtrlNet(curv)

