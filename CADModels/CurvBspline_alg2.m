%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the cubic B-spline curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% controlpoints
CtrlPts = zeros(2, 7);
CtrlPts(1 : 2, 1) = [0; 0];
CtrlPts(1 : 2, 2) = [0; 1];
CtrlPts(1 : 2, 3) = [1; 2];
CtrlPts(1 : 2, 4) = [2.5; -0.5];
CtrlPts(1 : 2, 5) = [4; 2];
CtrlPts(1 : 2, 6) = [5; 2.5];
CtrlPts(1 : 2, 7) = [6; 1];

% knot vector
KntVect = [0 0 0 0 1 2 3 4 4 4 4];
KntVect = KntVect ./ max(KntVect);

xi = linspace(0, 1, 101); % parametric points

% elvaluate the parametric points to plot the curve
C = BsplineEval({KntVect}, CtrlPts, {xi});

% elvaluate the knots
uqKntVect = KntVect([true; diff(KntVect(:)) > 0]); 
evalKnts = BsplineEval({KntVect}, CtrlPts, {uqKntVect});

figure
hold on
set(gcf, 'color', 'white');
axis off
axis equal
% Plot curve
plot(C(1, :), C(2, :), 'linewidth', 1.5);
% Plot control polygon
plot(CtrlPts(1, :), CtrlPts(2, :), 'k--');
% Plot control points
plot(CtrlPts(1, :), CtrlPts(2, :),...
    'r.','MarkerSize',15);
% Plot 5 points on the curve corresponding to 5 knots
plot(evalKnts(1, :), evalKnts(2, :),'s','MarkerEdgeColor', 'g',...
    'MarkerFaceColor', 'g', 'MarkerSize', 5);
axis equal
