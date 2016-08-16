%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a one-third circle modelled by a parameterized Bezier curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

radius = 1; % radius of the circle
alfa = 1/3 * pi; % 1/2 segment angle

% control points
CtrlPts = zeros(4, 3);
CtrlPts(1 : 2, 1) = [-sin(alfa); cos(alfa)];
CtrlPts(1 : 2, 2) = [0; 1/cos(alfa)];
CtrlPts(1 : 2, 3) = [sin(alfa); cos(alfa)];

CtrlPts(1 : 2, :) = CtrlPts(1 : 2, :) * radius;
% weights
CtrlPts(4, :) = 1;
fac = cos(alfa);
CtrlPts(:, 2) = CtrlPts(:, 2) * fac;

p = size(CtrlPts, 2) - 1; % degree of the curve

% evalutate the bernstein basis functions
xi = linspace(0, 1, 101); % allocate parametric points
B = AllBernstein(p, xi);

% contruct the rational bezier curve in 4D space
Cw = CtrlPts * B';

% project this curve to 3D space
w = Cw(4, :);
C = bsxfun(@rdivide, Cw, w);

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
% plot the curve
plot(C(1, :), C(2, :), 'b-', 'LineWidth', 1.5)
% plot control points
weights = CtrlPts(4, :);
plot(CtrlPts(1, :)./weights, CtrlPts(2, :)./weights, 'r.', 'MarkerSize', 20)
% plot control net
plot(CtrlPts(1, :)./weights, CtrlPts(2, :)./weights, 'k--')

axis off
axis equal





