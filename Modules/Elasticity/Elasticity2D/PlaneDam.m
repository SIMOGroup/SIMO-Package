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
H = 8;
D = 3;
alpha=pi/3;
L=D+H/tan(pi/3);

CtrlPts = zeros(4, 2, 2);

CtrlPts(1 : 3, 1, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2, 1) = [L; 0; 0];

CtrlPts(1 : 3, 1, 2) = [0; H; 0];
CtrlPts(1 : 3, 2, 2) = [D; H; 0];

CtrlPts(4, :, :) = 1;

KntVect{1} = [0 0 1 1];
KntVect{2} = [0 0 1 1];
% create NURBS surface in CAD 
Surf = CreateNURBS(KntVect, CtrlPts);
% degree of basis function
p=2;q=2;
% repeated knot value inside knot vector
kx=1;ky=1;
% number of elements in each direction
nelx=10; nely=10;
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
E = 105;
nu = 0.25;

KVals = calcLocalStiffnessMatrices2D(Mesh, Surf, E, nu, 'PlaneStrain');
[Rows, Cols, Vals] = convertToTripletStorage(Mesh, KVals);

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

f = zeros(Mesh.NDof, 1);
% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
P = 1;

Tx = @(x, y) P * (1 - y / H);

[Fx, DofsFx] = applyNewmannBdryVals(Surf, Mesh, Tx, 1, 'FX');

f(DofsFx) = f(DofsFx) + Fx;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
u_bar = @(x,y) 0;

[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, u_bar, 3, 'UX');
[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, u_bar, 3, 'UY');

BdryIdcs = [DofsY; DofsX];
BdryVals = [UY; UX];

FreeIdcs = setdiff(1 : Mesh.NDof, BdryIdcs);

d = zeros(Mesh.NDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

%------------------------------------------------------------------------
% strain energy
%------------------------------------------------------------------------
StrainEnergy = 0.5*f'*d;
disp(['StrainEnergy = ', num2str(StrainEnergy, 20)]);

% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 201), linspace(0, 1, 201)};
SPToVTKStress(Surf, d, ParaPts, E, nu, 'PlaneStrain', 'PlaneDamStress')
SPToVTK(Surf, d, ParaPts, 'PlaneDamDispl', 'Displ')

