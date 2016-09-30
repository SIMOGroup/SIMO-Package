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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points
tic;
a = 3;
b = 2;
CtrlPts = zeros(4, 2, 2);

CtrlPts(1 : 3, 1, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2, 1) = [a; 0; 0];

CtrlPts(1 : 3, 1, 2) = [0; b; 0];
CtrlPts(1 : 3, 2, 2) = [a; b; 0];

CtrlPts(4, :, :) = 1;

% knot vector

KntVect{1} = [0 0 1 1];
KntVect{2} = KntVect{1};

Surf = CreateNURBS(KntVect, CtrlPts);
Surf = KRefine(Surf, [5 4], [3 2], [0 0]);

figure
hold on
% grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
axis off
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'ScalarField');

% Assemble conduction matrix
disp([num2str(toc),'  Assembling the system'])
% material properties
ka = 1; % thermal conductivity
s = @(x, y) (6); %heat source
[KVals, FVals] = calcLocalConductionMatrices2D(Mesh, Surf, ka, s);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals);

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) sin(pi * x); 

[BCVals, BCIdcs] = projDrchltBdryVals(Surf, Mesh, h, 4, 'TEMP');

FreeIdcs = 1 : Mesh.NDof;
HomoBCIdcs = union(Mesh.Boundary(1).Dofs, Mesh.Boundary(2).Dofs);
HomoBCIdcs = union(HomoBCIdcs, Mesh.Boundary(3).Dofs);
FreeIdcs(union(HomoBCIdcs, BCIdcs)) = [];

d = zeros(Mesh.NDof, 1);
d(BCIdcs) = BCVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BCIdcs) * BCVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

% % % Plot heat flux at gauss points
% % disp([num2str(toc),'  Ploting heat flux'])
% % PlotFlux2D(Surf, Mesh, d, ka)

% Evaluate and plot analytical solution vs approximate solution
disp([num2str(toc),'  Evaluating analytical solution'])
[C, F] = NURBSEval(Surf, {linspace(0, 1, 101), 1}, d);
x = C(1, :);
y = C(2, :);

fx = @(x) sin(pi * x);
laplace = SPCondRectExactLaplace(fx, x, y, a, b, 20);
poisson = SPCondRectExactPoisson(s, x, y, a, b, 20, 20);
TExact = laplace + poisson;

figure
hold on
plot(x, TExact, 'r', 'LineWidth', 1.5);
plot(x, F, 'b');
legend('exact', 'approx')

% Export to *.vts filez
disp([num2str(toc),'  Export result to *.vts file'])
ParaPts = {linspace(0, 1, 31), linspace(0, 1, 31)};
SPToVTK(Surf, d, ParaPts, 'SPCondRect', 'Temperature')
clear sym


