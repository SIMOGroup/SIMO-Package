clear
close all
clc

KntVect = [0 0 0 1 2 3 4 5 6 7 8 9 10 10 10];
KntVect = KntVect ./ max(KntVect);
p = 2; % order of basis functions
NCtrlPts = numel(KntVect) - p - 1;

xi = linspace(0, 1, 1001);
Idx = FindSpan(NCtrlPts, p, xi, KntVect);

N = DersBasisFuns(Idx, xi, p, 1, KntVect);
N0 = squeeze(N(:, :, 1));

figure 
hold on
set(gcf, 'color', 'white');
plot(xi, N0');




