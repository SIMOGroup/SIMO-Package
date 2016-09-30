% For testing purpose

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
format short

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points
tic
r = 0.02;
t = 0.08;
h = 0.14;

Curv{1} = NURBSCirc(r, [0, 0, 0], 0, pi/2);
Curv{2} = NURBSCirc(r + t, [0, 0, 0], 0, pi/2);

Surf = NURBSRuled(Curv{1}, Curv{2});

Volu = NURBSExtrude(Surf, [0 0 h]);

Volu = KRefine(Volu, [8 8 10], [2 4 4], [1 1 1]);

figure
hold on
% grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
% axis off
view(3)
PlotGeo(Volu);
PlotKnts(Volu);
PlotCtrlPts(Volu);
PlotCtrlNet(Volu);

% --------------------------------------------------------------
Mesh3d = Mesh3D(Volu, 'ScalarField');

% Assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])
% material properties
k = 52; % thermal conductivity
s = @(x, y, z) (0); %heat source

[KVals, FVals] = calcLocalConductionMatrices3D(Mesh3d, Volu, k, s);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh3d, KVals, FVals);

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
g = @(x, y, z) ((z - 0.04) >= 1e-6 && (0.1 - z) >= 1e-6) * 5e5;
Refs = 3;
[flux, Dofs] = applyNewmannBdryVals(Volu, Mesh3d, g, Refs, 'HFLUX');

f(Dofs) = f(Dofs) + flux;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
T = 0; % degree Celsius
% T = 273.15; % degree Kelvin

BdryIdcs = unique([Mesh3d.Boundary(4).Dofs, Mesh3d.Boundary(5).Dofs, Mesh3d.Boundary(6).Dofs]);
BdryVals = repmat(T, numel(BdryIdcs), 1);

FreeIdcs = setdiff(1 : Mesh3d.NDof, BdryIdcs);

d = zeros(Mesh3d.NDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])

ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101), linspace(0, 1, 201)};
SPToVTK(Volu, d, ParaPts, 'SPCondCylind', 'Temperature')

[C, benchmark] = NURBSEval(Volu, {0, 0.02/0.08, 0.04/0.14}, d)

