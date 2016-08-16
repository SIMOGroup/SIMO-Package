close all
clear
clc

% control points

CtrlPts = zeros(4, 4);
CtrlPts(1 : 3, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2) = [2; 0; 0];
CtrlPts(1 : 3, 3) = [1; sqrt(3); 0];
CtrlPts(1 : 3, 4) = [3; sqrt(3); 0];

CtrlPts(4, :) = 1;

W = [0.1, 0.5, 1, 3];
Knts = [0 0 0 0.5 1 1 1];

Curv = cell(numel(W), 1);

for i = 1 : numel(W)   
    temp = CtrlPts;
    temp(:, 2) = temp(:, 2) .* W(i);
    Curv{i} = CreateNURBS({Knts}, temp);
end

figure 
hold on
axis off
set(gcf,'color','white')
set(gca,'XTick', 0 : 1: 3)
set(gca,'YTick', 0 : 1: 2)

for i = 1 : numel(W)
    PlotGeo(Curv{i})
    PlotCtrlPts(Curv{i})
    PlotKnts(Curv{i})
    PlotCtrlNet(Curv{i})
end

