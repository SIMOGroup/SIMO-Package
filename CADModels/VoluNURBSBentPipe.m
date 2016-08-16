clc
clear
close all

Curv{1} = NURBSCirc(1);
Curv{2} = NURBSCirc(2);

Annulus = NURBSRuled(Curv{1}, Curv{2});

Pipe = NURBSExtrude(Annulus, [0 0 3]);
Elbow = NURBSRevolve(Annulus, [3, 0, 0], [0, -1, 0], pi / 2);

Elbow = NURBSPermute(Elbow, [2 3 1]);
Pipe = NURBSReverse(Pipe, 3);

bentPipe = NURBSJoin(Pipe, Elbow, 3);

figure
hold on
grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
view(3)
PlotGeo(bentPipe)
PlotKnts(bentPipe)
PlotCtrlPts(bentPipe)
PlotCtrlNet(bentPipe)
