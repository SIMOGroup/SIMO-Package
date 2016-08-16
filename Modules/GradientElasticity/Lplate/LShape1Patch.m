clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
% control points
tic;
L = 1;
D = 0.5;

% patch 1
CtrlPts = zeros(4, 3, 2);

CtrlPts(1 : 3, 1, 1) = [0; L; 0];
CtrlPts(1 : 3, 2, 1) = [0; 0; 0];
CtrlPts(1 : 3, 3, 1) = [L; 0; 0];

CtrlPts(1 : 3, 1, 2) = [D; L; 0];
CtrlPts(1 : 3, 2, 2) = [D; D; 0];
CtrlPts(1 : 3, 3, 2) = [L; D; 0];

CtrlPts(4, :, :) = 1;

% knot vector

KntVect1 = [0 0 0.5 1 1];
KntVect2 = [0 0 1 1];

Surf = CreateNURBS({KntVect1, KntVect2}, CtrlPts);

figure
hold on
axis equal
daspect([1, 1, 1])
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf);
PlotCtrlNet(Surf);


