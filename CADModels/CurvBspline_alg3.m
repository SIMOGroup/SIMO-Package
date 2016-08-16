clc
close all
clear

% control points
CtrlPts = zeros(4, 7);
CtrlPts(1 : 2, 1) = [1; 0];
CtrlPts(1 : 2, 2) = [0.5; 1.5];
CtrlPts(1 : 2, 3) = [2; 2.5];
CtrlPts(1 : 2, 4) = [3; 1];
CtrlPts(1 : 2, 5) = [2.5; 0];
CtrlPts(1 : 2, 6) = [4.5; 0.5];
CtrlPts(1 : 2, 7) = [3.5; 1.5];

CtrlPts(4, :) = 1;

KntVect{1} = [0 0 0 0 0 1 2 3 3 3 3 3];
KntVect{1} = KntVect{1} ./ max(KntVect{1});

Curv = CreateNURBS(KntVect, CtrlPts);

figure
hold on
axis off
set(gcf,'color','white')
set(gca,'XTick', 0 : 1: 3)
set(gca,'YTick', 0 : 1: 2)

PlotGeo(Curv)
PlotCtrlPts(Curv)
PlotKnts(Curv)
PlotCtrlNet(Curv)
