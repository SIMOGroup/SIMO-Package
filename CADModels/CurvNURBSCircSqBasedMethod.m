%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% representing a full circle using nurbs curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

radius = 1;
% control points for a circle
CtrlPts = zeros(4, 9);

CtrlPts(1 : 2, 1, 1) = [1; 0];
CtrlPts(1 : 2, 2, 1) = [1; 1];
CtrlPts(1 : 2, 3, 1) = [0; 1];
CtrlPts(1 : 2, 4, 1) = [-1; 1];
CtrlPts(1 : 2, 5, 1) = [-1; 0];
CtrlPts(1 : 2, 6, 1) = [-1; -1];
CtrlPts(1 : 2, 7, 1) = [0; -1];
CtrlPts(1 : 2, 8, 1) = [1; -1];
CtrlPts(1 : 2, 9, 1) = [1; 0];

% scale radius of circle
CtrlPts(1 : 3, :, :) = CtrlPts(1 : 3, :, :) * radius;

% assign the weights
CtrlPts(4, :, :) = 1;
W = 1/sqrt(2);
%W = 1 / sqrt(2); % W_i = 1/sqrt(2) if i even
CtrlPts(:, 2 : 2 : end, :) = CtrlPts(:, 2 : 2 : end, :) * W;

% knot vector
KntVect = [0 0 0 1 1 2 2 3 3 4 4 4];
KntVect = KntVect ./ max(KntVect);

ParaPts = linspace(0, 1, 101); % parametric points

NCtrlPts = size(CtrlPts, 2);
p = numel(KntVect) - NCtrlPts - 1;

Idx = FindSpan(NCtrlPts, p, ParaPts, KntVect);
N0 = BasisFuns(Idx, ParaPts, p, KntVect);

Cw = zeros(4, numel(ParaPts));
for i = 1 : p + 1
    Cw = Cw + bsxfun(@times, N0(:, i)', CtrlPts(:, Idx - p + i - 1));
end

% project the curve into Cartesian 3D space
C = Cw(1 : 3, :);
W = Cw(4, :);
C = bsxfun(@rdivide, C, W);

% project the control points into Cartesian 3D space
CtrlPts3D = bsxfun(@rdivide, CtrlPts(1 : 3, :), CtrlPts(4, :));

figure
hold on
axis off
set(gcf, 'color', 'white');
daspect([1 1 1]);
axis equal
% Plot curve
plot(C(1, :), C(2, :), 'linewidth', 1.5);
% Plot control polygon
plot(CtrlPts3D(1, :), CtrlPts3D(2, :), 'k--');
% Plot control points
plot(CtrlPts3D(1, :), CtrlPts3D(2, :), 'r.','MarkerSize',20);
