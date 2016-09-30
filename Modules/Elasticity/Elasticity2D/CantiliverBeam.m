% CantileverBeam
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

close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points
tic
L = 48;
D = 12;
CtrlPts = zeros(4, 2, 2);

CtrlPts(1 : 3, 1, 1) = [0; -D/2; 0];
CtrlPts(1 : 3, 2, 1) = [L; -D/2; 0];

CtrlPts(1 : 3, 1, 2) = [0; D/2; 0];
CtrlPts(1 : 3, 2, 2) = [L; D/2; 0];

CtrlPts(4, :, :) = 1;

KntVect{1} = [0 0 1 1];
KntVect{2} = [0 0 1 1];
% create NURBS surface in CAD 
Surf = CreateNURBS(KntVect, CtrlPts);
% degree of basis function
p=4;q=2;
% repeated knot value inside knot vector
kx=1;ky=1;
% number of elements in each direction
nelx=10; nely=2;
% h,p,k-refinements
Surf = KRefine(Surf, [nelx, nely], [p, q], [p-kx, q-ky]);

figure
hold on
axis equal
daspect([1, 1, 1])
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf);
PlotCtrlNet(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'VectorField');

% Assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])
% material properties
E = 3e7;
nu = 0.3;

KVals = calcLocalStiffnessMatrices2D(Mesh, Surf, E, nu, 'PlaneStress');
[Rows, Cols, Vals] = convertToTripletStorage(Mesh, KVals);

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

f = zeros(Mesh.NDof, 1);
% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
P = 1000;
I = D ^ 3 / 12;

Tauy = @(x, y) -P * (D ^ 2 / 4 - y ^ 2) / (2 * I);

[Fy, DofsFy] = applyNewmannBdryVals(Surf, Mesh, Tauy, 2, 'FY');

f(DofsFy) = f(DofsFy) + Fy;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
solu = CantiliverBeamExactSolu(L, D, E, I, nu, P);

[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, solu.ux, 1, 'UX');
[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, solu.uy, 1, 'UY');

BdryIdcs = [DofsY; DofsX];
BdryVals = [UY; UX];

FreeIdcs = setdiff(1 : Mesh.NDof, BdryIdcs);

d = zeros(Mesh.NDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);
% break
%------------------------------------------------------------------------
% strain energy
%------------------------------------------------------------------------
StrainEnergy = 0.5*f'*d;
disp(['StrainEnergy = ', num2str(StrainEnergy, 20)]);
% compare stress between IGA and exact solution
% -----------------------------------------------------------------------
% IGA solution
[C, F] = NURBSEval(Surf, {linspace(0, 1, 101), 0.5}, d);

x = C(1, :)';
y = C(2, :)';

figure
hold on
title('Vertical displacement at the center line')
plot(x, solu.uy(x, y), 'r-');
plot(x, F(2, :)', 'b*');
legend('exact', 'approx');
% -----------------------------------------------------------------------

% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 201), linspace(0, 1, 201)};
SPToVTKStress(Surf, d, ParaPts, E, nu, 'PlaneStress', 'CantiliverBeamStress')
SPToVTK(Surf, d, ParaPts, 'CantiliverBeamDispl', 'Displ')

