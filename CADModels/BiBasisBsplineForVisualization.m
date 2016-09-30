% -------------------------------------------------------------------------
% Represent bivariate B-spline basis functions using "OneBasisFun" routine
% (for visualization purpose)
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

KntVect{1} = [0 0 0 0 1 2 3 4 4 4 4];
KntVect{1} = KntVect{1} ./ max(KntVect{1});

KntVect{2} = [0 0 0 1 2 3 3 4 5 5 5];
KntVect{2} = KntVect{2} ./ max(KntVect{2});

p = 3; % order of basis functions in xi direction
q = 2; % order of basis functions in eta direction

[X, Y] = meshgrid(linspace(0, 1, 41));
R1 = zeros(size(X, 1), size(X, 1));
R2 = zeros(size(X, 1), size(X, 1));

for i = 1 : size(X, 1)
    for j = 1 : size(X, 1)
        eta = X(1, i);
        xi = Y(j, 1);
        N = OneBasisFun(p, KntVect{1}, 5, xi);
        M = OneBasisFun(q, KntVect{2}, 3, eta);
        R1(i, j) = N * M;
    end
end

for i = 1 : size(X, 1)
    for j = 1 : size(X, 1)
        eta = X(1, i);
        xi = Y(j, 1);
        N = OneBasisFun(p, KntVect{1}, 3, xi);
        M = OneBasisFun(q, KntVect{2}, 8, eta);
        R2(i, j) = N * M;
    end
end

% Data for plot univariate basis functions
ParaPts = 201;
xi = linspace(0, 1, ParaPts);

nen(1) = numel(KntVect{1}) - p - 1; % number of local basis function
N = zeros(numel(xi), nen(1));
for i = 1 : nen(1)
    N(:, i) = OneBasisFun(p, KntVect{1}, i, xi);
end

eta = linspace(0, 1, ParaPts);
nen(2) = numel(KntVect{2}) - q - 1; % number of local basis function
M = zeros(numel(eta), nen(2));
for i = 1 : nen(2)
    M(:, i) = OneBasisFun(q, KntVect{2}, i, eta);
end

% Plot bivariate basis function
figure
hold on
set(gcf,'color','white')
set(gca,'XTick', 0:1/4:1)
set(gca,'YTick', 0:1/5:1)
set(gca,'ZTick', 0:0.25:1)
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
axis off
axis equal
R = R1 + R2;
surf(X, Y, R)
colormap(autumn)
xlabel('x')
ylabel('y')
zlabel('z')
plot3(xi, zeros(ParaPts, 1), N);
plot3(zeros(ParaPts, 1), eta, M);
view(3)
