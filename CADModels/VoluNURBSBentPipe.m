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
