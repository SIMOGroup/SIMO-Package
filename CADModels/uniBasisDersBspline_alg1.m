%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Bspline basis functions and their first derivarives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
p = 2; % order of basis functions
KntVect = [zeros(1,p+1) 1/2 ones(1,p+1)];
%KntVect = [0 0 1 1]; % knot vector
%KntVect = [0 0 0 1 2 3 4 4 4]; % knot vector
KntVect = KntVect ./ max(KntVect);
%p = 2; % order of basis functions
NCtrlPts = numel(KntVect) - p - 1; % number of control points

xi = linspace(0, 1, 1000);

N = zeros(NCtrlPts, numel(xi), p + 1);
for i = 1 : NCtrlPts
    N(i, :, :) = DersOneBasisFun(p, KntVect, i, xi, p);
end

% Bspline basis functions
N0 = squeeze(N(:, :, 1));

% First derivatives of Bspline basis functions
N1 = squeeze(N(:, :, 2));

figure
set(gcf, 'color', 'white');
plot(xi, N0);

figure
set(gcf, 'color', 'white');
plot(xi, N1);


