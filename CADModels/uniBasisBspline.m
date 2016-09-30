% -------------------------------------------------------------------------
% Plot B-spline basis functions
% -------------------------------------------------------------------------

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

