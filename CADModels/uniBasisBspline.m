%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Bspline basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

p = 2; % order of basis functions
KntVect = [zeros(1,p+1) ones(1,p+1)];
%KntVect = [0 0 0 1 2 3 4 5 6 7 8 9 10 10 10];
KntVect = KntVect ./ max(KntVect);
NCtrlPts = numel(KntVect) - p - 1; % number of control points
ParaPts = linspace(0, 1, 1000); % parametric points \xi in parameter space

% % % Separate basis functions
% % % i.e. N0(:, i) is the ith basis function
% % % size(N0) = [numel(ParaPts), NCtrlPts]
% % N0 = zeros(numel(ParaPts), NCtrlPts);
% % for i = 1 : NCtrlPts
% %     N0(:, i) = OneBasisFun(p, KntVect, i, ParaPts);
% % end

Idx = FindSpan(NCtrlPts, p, ParaPts, KntVect);
N0 = BasisFuns(Idx, ParaPts, p, KntVect);

figure
set(gcf, 'color', 'white');
% plot(ParaPts, N0);
plot(ParaPts, N0(:, 1 : p + 1));

