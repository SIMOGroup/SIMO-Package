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
