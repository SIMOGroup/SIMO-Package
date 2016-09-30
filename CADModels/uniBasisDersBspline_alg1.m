% -------------------------------------------------------------------------
% Plot Bspline basis functions and their first derivarives
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


