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
R1 = 2;
R2 = 0.75;

Curv{1} = NURBSCirc(R2);
Curv{2} = NURBSCirc(R1);

Surf = NURBSRuled(Curv{1}, Curv{2});

Surf = KRefine(Surf, [10 4], [5 5], [4 4]);
% Surf = KRefine(Surf, [2 2], [2 3], [1 2]);

% Material properties

E  = 3e6;
nu = 0.3;
t  = 0.01; % thickness

% Boundary condition

q = -1;

figure
hold on
% grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
axis off
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf, 1);
% PlotCtrlNet(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'Plate');
% material properties

% Coupling coincide control points by global numbering
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
[KVals, FVals] = calcLocalStiffnessMatricesKirchhoffPlate(Mesh, Surf, t, E, nu, q);
[Rows, Cols, Vals, ValsF] = convertToTripletStorage(Mesh, KVals, FVals);

Rs = GNum(Rows);
Cs = GNum(Cols);
% Convert triplet data to sparse matrix
K = sparse(Rs, Cs, Vals);

f = accumarray(GNum(1 : numel(ValsF)), ValsF);

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;

cp1 = Mesh.Boundary(1).Dofs;
cp2 = Mesh.Boundary(1).NextLayerDofs;
cp3 = Mesh.Boundary(2).NextLayerDofs;
coupledNodes1 = [cp3 cp1; cp1 cp2];
coupledNodes2 = findConn(Surf, 0.25, 1);
coupledNodes3 = findConn(Surf, 0.5, 1);
coupledNodes4 = findConn(Surf, 0.75, 1);

w = 1e7;
penaltyStiffness = w*[1 -1;-1 1];
coupledNodes = [coupledNodes1; coupledNodes2; coupledNodes3; coupledNodes4];
for i = 1 : size(coupledNodes, 1)
    sctr  = GNum(coupledNodes(i, :))';

    K(sctr,sctr) = K(sctr,sctr) + penaltyStiffness;
end

% Inner curve
[u3, Dofs3] = projDrchltBdryVals(Surf, Mesh, h, 3, 'PLATE', GNum);
% Outer curve
[u4, Dofs4] = projDrchltBdryVals(Surf, Mesh, h, 4, 'PLATE', GNum);

NextDofs3 = GNum(Mesh.Boundary(3).NextLayerDofs);
NextDofs4 = GNum(Mesh.Boundary(4).NextLayerDofs);

BCIdx = unique([Dofs3; NextDofs3; Dofs4; NextDofs4]);
BCVal = [u3; zeros(numel(u3), 1); u4; zeros(numel(u4), 1)];

FreeIdx = setdiff(1 : GDof, BCIdx);

d = zeros(GDof, 1);
d(BCIdx) = BCVal;

f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdx) * BCVal; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdx) = K(FreeIdx, FreeIdx) \ f(FreeIdx);

d = d(GNum);
dd = [zeros(2 * numel(d), 1); d];

% Export result to *.vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])

ParaPts = {linspace(0, 1, 101), linspace(0, 1, 51)};
C = NURBSEval(Surf, ParaPts);
x = squeeze(C(1, :, :));
y = squeeze(C(2, :, :));

solu = SPAnnularPlateExactSolu(R2, R1, t, E, nu, q);
uz = solu.w(x, y);
duz = solu.theta(x, y);
uuz = zeros(size(C));
uuz(3, :, :) = uz;


exportToVTK(C, squeeze(uuz), 'SPAnnularPlateExact', 'Displ');
SPToVTK(Surf, dd, ParaPts, 'SPAnnularPlate', 'Displ')

