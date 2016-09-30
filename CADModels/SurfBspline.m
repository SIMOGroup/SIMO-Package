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

close all
clear
clc

% control points

CtrlPts = zeros(3, 4, 5);

CtrlPts(1 : 3, 1, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2, 1) = [1/3; 0; 0];
CtrlPts(1 : 3, 3, 1) = [2/3; 0; 0];
CtrlPts(1 : 3, 4, 1) = [1; 0; 0];

CtrlPts(1 : 3, 1, 2) = [0; 1/4; 0];
CtrlPts(1 : 3, 2, 2) = [1/3; 1/4; 1/4];
CtrlPts(1 : 3, 3, 2) = [2/3; 1/4; 0];
CtrlPts(1 : 3, 4, 2) = [1; 1/4; 0];

CtrlPts(1 : 3, 1, 3) = [0; 2/4; 0];
CtrlPts(1 : 3, 2, 3) = [1/3; 2/4; 1/4];
CtrlPts(1 : 3, 3, 3) = [2/3; 2/4; 0];
CtrlPts(1 : 3, 4, 3) = [1; 2/4; 0];

CtrlPts(1 : 3, 1, 4) = [0; 3/4; 0];
CtrlPts(1 : 3, 2, 4) = [1/3; 3/4; 0];
CtrlPts(1 : 3, 3, 4) = [2/3; 3/4; 0];
CtrlPts(1 : 3, 4, 4) = [1; 3/4; 0];

CtrlPts(1 : 3, 1, 5) = [0; 1; 0];
CtrlPts(1 : 3, 2, 5) = [1/3; 1; 0];
CtrlPts(1 : 3, 3, 5) = [2/3; 1; 0];
CtrlPts(1 : 3, 4, 5) = [1; 1; 0];

% knot vectors
KntVect{1} = [0 0 0 1/2 1 1 1];
KntVect{2} = [0 0 0 1/3 2/3 1 1 1];

figure
hold on
set(gcf,'color','white')
xlabel('x');
ylabel('y');
zlabel('z');
daspect([1 1 1])
axis equal
axis tight
axis off
rotate3d on

% Plot the control net
% Plot control points along y-direction
for i = 1 : size(CtrlPts, 2)
    Pts = reshape(CtrlPts(1 : 3, i, :), 3, []);
    plot3(Pts(1, :), Pts(2, :), Pts(3, :),'k--')
    plot3(Pts(1, :), Pts(2, :), Pts(3, :),'r.','MarkerSize', 15)
end
% Plot control points along x-direction
for j = 1 : size(CtrlPts, 3)
    Pts = reshape(CtrlPts(1 : 3, :, j), 3, []);
    plot3(Pts(1, :), Pts(2, :), Pts(3, :),'k--')
    plot3(Pts(1, :), Pts(2, :), Pts(3, :),'r.','MarkerSize', 15)
end 
% Plot the B-spline surface

ParaPts{1} = linspace(0, 1, 101);
ParaPts{2} = linspace(0, 1, 101);

S = BsplineEval(KntVect, CtrlPts, ParaPts);

bcol = [173 234 234] ./ 250;
bcol = repmat(bcol, 3, 1);
colormap(bcol);
surfl(squeeze(S(1,:,:)), squeeze(S(2,:,:)), squeeze(S(3,:,:)));
shading interp
view(3)

% Plot the knots
uqKntVect{1} = unique(KntVect{1});
uqKntVect{2} = unique(KntVect{2});
Curv{1} = BsplineEval(KntVect, CtrlPts, {uqKntVect{1}, ParaPts{2}});
Curv{2} = BsplineEval(KntVect, CtrlPts, {ParaPts{1}, uqKntVect{2}});

for i = 1 : numel(uqKntVect{1})
    EtaCurv = reshape(Curv{1}(:, i, :), 3, []);
    plot3(EtaCurv(1, :), EtaCurv(2, :), EtaCurv(3, :), 'k');
end
for i = 1 : numel(uqKntVect{2})
    XiCurv = reshape(Curv{2}(:, :, i), 3, []);
    plot3(XiCurv(1, :), XiCurv(2, :), XiCurv(3, :), 'k');
end
