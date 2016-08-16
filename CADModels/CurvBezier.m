%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An example of a Bezier curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

% control points
CtrlPts = zeros(2, 3);
CtrlPts(:, 1) = [1; -1];
CtrlPts(:, 2) = [2; 1];
CtrlPts(:, 3) = [3; 0];

p = size(CtrlPts, 2) - 1; % degree of the curve

% evalutate the bernstein basis functions
xi = linspace(0, 1, 101);
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

% plot the curve, control points and control polygon
figure
hold on
set(gcf, 'color', 'white')
plot(C(:,1), C(:,2), 'b-', 'LineWidth', 1.5)
plot(CtrlPts(1, :), CtrlPts(2, :), 'r.', 'MarkerSize', 15)
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--')
axis off
axis equal
