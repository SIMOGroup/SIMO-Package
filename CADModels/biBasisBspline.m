% -------------------------------------------------------------------------
% Represent bivariate B-spline basis functions
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

Xi = [0 0 0 0 1 2 3 4 4 4 4];
Xi = Xi ./ max(Xi);

Eta = [0 0 0 1 2 3 4 5 5 5];
Eta = Eta ./ max(Eta);

p = 3; % order of basis functions in xi direction
q = 2; % order of basis functions in eta direction

NCtrlPtsXi = numel(Xi) - p - 1;
NCtrlPtsEta = numel(Eta) - q - 1;

NPts = [401 201];
x = linspace(0, 1, NPts(1));
y = linspace(0, 1, NPts(2));

IdxXi = FindSpan(NCtrlPtsXi, p, linspace(0, 1, NPts(1)), Xi);
N0Xi = BasisFuns(IdxXi, linspace(0, 1, NPts(1)), p, Xi);

IdxEta = FindSpan(NCtrlPtsEta, q, linspace(0, 1, NPts(2)), Eta);
N0Eta = BasisFuns(IdxEta, linspace(0, 1, NPts(2)), q, Eta);

N0 = zeros(NPts(1), NPts(2), (p + 1) * (q + 1));
for i = 1 : NPts(1)
    for j = 1 : NPts(2)
        tmp = N0Xi(i, :)' * N0Eta(j, :);
        N0(i, j, :) = tmp(:);
    end
end

% Plot bivariate basis function
figure
hold on
set(gcf,'color','white')
axis off
axis equal
for i = 1 : (p + 1) * (q + 1)
    surfl(repmat(x', 1, NPts(2)), repmat(y, NPts(1), 1), N0(:, :, i))
end
%shading interp
%colormap(summer)
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
plot3(x, zeros(NPts(1), 1), N0Xi);
plot3(zeros(NPts(2), 1), y, N0Eta);
