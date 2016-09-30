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
R1 = 4;
R2 = 2;

Curv{1} = NURBSCirc(R2);
Curv{2} = NURBSCirc(R1);

Surf = NURBSRuled(Curv{1}, Curv{2});

% Surf = KRefine(Surf, [2 1], [2 1], [0 0]);
Surf = KRefine(Surf, [10 10], [4 4], [0 0]);

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
PlotCtrlNet(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'ScalarField');
% material properties
ka = 1; % thermal conductivity - W/(m * ^\circ C)
s = @(x, y) (-(12 * y ^ 2 + 6 * x));

% Coupling coincident control points by global numbering
GNum = zeros(Mesh.NDof, 1);
gluedDofs = union(Mesh.Boundary(1).Dofs, Mesh.Boundary(2).Dofs);
nonGludedDofs = setdiff(1 : Mesh.NDof, gluedDofs);
GNum(nonGludedDofs) = 1 : numel(nonGludedDofs);

newDofs = numel(nonGludedDofs) + (1 : numel(Mesh.Boundary(1).Dofs));
GNum(Mesh.Boundary(1).Dofs) = newDofs;
GNum(Mesh.Boundary(2).Dofs) = newDofs;

GDof = numel(nonGludedDofs) + numel(newDofs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------PROCESSING--------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'  Assembling the system'])
[KVals, FVals] = calcLocalConductionMatrices2D(Mesh, Surf, ka, s);
[Rows, Cols, Vals, ValsF] = convertToTripletStorage(Mesh, KVals, FVals);

Rs = GNum(Rows);
Cs = GNum(Cols);
% Convert triplet data to sparse matrix
K = sparse(Rs, Cs, Vals);

f = accumarray(GNum(1 : numel(ValsF)), ValsF);

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) x ^ 3 + y ^ 4;

% Inner curve
[BCVals, BCIdcs] = projDrchltBdryVals(Surf, Mesh, h, [3, 4], 'TEMP', GNum);

FreeIdx = setdiff(1 : GDof, BCIdcs);

d = zeros(GDof, 1);
d(BCIdcs) = BCVals;

f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdcs) * BCVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdx) = K(FreeIdx, FreeIdx) \ f(FreeIdx);

d = d(GNum);

% % % Plot heat flux at gauss points
% % disp([num2str(toc),'  Ploting heat flux'])
% % PlotFlux2D(Surf, Mesh, d, ka)

[C, F] = NURBSEval(Surf, {linspace(0, 0.25, 101), 0.9}, d);

alfa = atan2(C(2, :), C(1, :));
plot(C(1, :), C(2, :));

figure
hold on
x = C(1, :);
y = C(2, :);
plot(rad2deg(alfa), x.^3 + y.^4, 'r', 'LineWidth', 1.5);
plot(rad2deg(alfa), F, 'b')

% Export result to *.vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])

ParaPts = {linspace(0, 1, 101), linspace(0, 1, 51)};
SPToVTK(Surf, d, ParaPts, 'SPCondAnnular', 'Temperature')

CC = NURBSEval(Surf, ParaPts);
x = CC(1, :, :);
y = CC(2, :, :);
exact = x.^3 + y.^4;
exportToVTK(CC, squeeze(exact), 'SPCondAnnularExact', 'Temperature');
