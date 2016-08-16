
clear
close all
clc

% control points in homogeneous 4D coordinates
CtrlPts = zeros(2, 7);
CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [0; 1];
CtrlPts(1 : 2, 3) = [1; 2];
CtrlPts(1 : 2, 4) = [2.5; -0.5];
CtrlPts(1 : 2, 5) = [4; 2];
CtrlPts(1 : 2, 6) = [5; 2.5];
CtrlPts(1 : 2, 7) = [6; 1];

p = size(CtrlPts, 2) - 1; % degree of the curve

% evalutate the bernstein basis functions
xi = linspace(0, 1, 401);
B = AllBernstein(p, xi);

% contruct bezier curve
C = B * CtrlPts';

% plot basis functions
figure
set(gcf,'color','white')
axis equal
daspect([1 1 1])
set(gcf, 'color', 'white');
plot(xi, B, 'LineWidth', 1.5)

% plot curve, control points and control polygon
figure
hold on
set(gcf, 'color', 'white')
% project to cartesian 3D coordinates
plot(C(:,1), C(:, 2), 'b-', 'LineWidth', 1.5)
plot(CtrlPts(1, :), CtrlPts(2, :),'r.', 'MarkerSize', 15)
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--')
axis equal
axis off

