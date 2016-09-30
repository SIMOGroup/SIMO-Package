% Sheet with Small Circular Hole

close all
clear
clc
addpath('C:\Program Files\Mosek\7\toolbox\r2013aom');

%{
Copyright (C) <2014-2016>  <Hung Nguyen-Xuan, Khanh Chau-Nguyen>

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

tic
% define problem 
% R = 0.25; % m
% L = 4 * R;

% control points
CtrlPts = zeros(4, 6, 5);
CtrlPts(1 : 3, 1, 1) = [1; 0; 0];
CtrlPts(1 : 3, 2, 1) = [1.5; 0; 0];
CtrlPts(1 : 3, 3, 1) = [1.75; 0; 0];
CtrlPts(1 : 3, 4, 1) = [2.25; 0; 0];
CtrlPts(1 : 3, 5, 1) = [2.5; 0; 0];
CtrlPts(1 : 3, 6, 1) = [3.0; 0; 0];

CtrlPts(1 : 3, 1, 2) = [0.99; 0.26; 0];
CtrlPts(1 : 3, 2, 2) = [1.44; 0.77; 0];
CtrlPts(1 : 3, 3, 2) = [1.7; 0.74; 0];
CtrlPts(1 : 3, 4, 2) = [2.25; 0.75; 0];
CtrlPts(1 : 3, 5, 2) = [2.55; 0.78; 0];
CtrlPts(1 : 3, 6, 2) = [3.0; 0.26; 0];

CtrlPts(1 : 3, 1, 3) = [0.73; 0.73; 0];
CtrlPts(1 : 3, 2, 3) = [0.91; 2.32; 0];
CtrlPts(1 : 3, 3, 3) = [1.5; 2.22; 0];
CtrlPts(1 : 3, 4, 3) = [2.5; 2.2; 0];
CtrlPts(1 : 3, 5, 3) = [3.07; 2.31; 0];
CtrlPts(1 : 3, 6, 3) = [3.26; 0.73; 0];

CtrlPts(1 : 3, 1, 4) = [0.26; 0.99; 0];
CtrlPts(1 : 3, 2, 4) = [0.32; 3.49; 0];
CtrlPts(1 : 3, 3, 4) = [0.74; 3.45; 0];
CtrlPts(1 : 3, 4, 4) = [3.25; 3.46; 0];
CtrlPts(1 : 3, 5, 4) = [3.67; 3.49; 0];
CtrlPts(1 : 3, 6, 4) = [3.73; 0.99; 0];

CtrlPts(1 : 3, 1, 5) = [0; 1; 0];
CtrlPts(1 : 3, 2, 5) = [0; 4; 0];
CtrlPts(1 : 3, 3, 5) = [0; 4; 0];
CtrlPts(1 : 3, 4, 5) = [4; 4; 0];
CtrlPts(1 : 3, 5, 5) = [4; 4; 0];
CtrlPts(1 : 3, 6, 5) = [4; 1; 0];

CtrlPts(4, :) = 1;

KntVect{1} = [0 0 0 1/4 2/4 3/4 1 1 1];
KntVect{2} = [0 0 0 1/3 2/3 1 1 1];

Surf = CreateNURBS(KntVect, CtrlPts);

% degree of basis function
p=3;q=3;
% repeated knot value inside knot vector
kx=1;ky=1;
% number of elements in each direction
nelx=4; nely=4;
% h,p,k-refinements
Surf = KRefine(Surf, [nelx, nely], [p, q], [p-kx, q-ky]);

% show geometry
figure
hold on
axis equal
daspect([1, 1, 1])
axis off
PlotGeo(Surf);
PlotKnts(Surf);
%PlotCtrlPts(Surf);
%PlotCtrlPts(Surf, 1);
%PlotCtrlNet(Surf);
% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'VectorField');

% assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])

% material properties
sigmap = 1; % Yield stress MPa
P1 = 1;     % load

%stressState ='PlaneStress';
stressState ='PlaneStrain';

% define problem type and total number of DOFs of system
if ( strcmp(stressState, 'PlaneStress') )       % Plane stress case
    noDofs=Mesh.NDof;
else
    imode=0;
    noDofs=Mesh.NDof+2*Mesh.NEl*imode;
end

%[Bstrain, WeJ, countGP] = calcLocalBMatrices2D(Mesh, Surf, stressState);

[Bstrain, WeJ, countGP] = calcLocalBMatrices2DOp(Mesh, Surf, stressState);

Wex = zeros(noDofs, 1);

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])

Ty = @(x, y) P1 * (abs(y - 4) <= 1e-6);

[Fy, DofsFy] = applyNewmannBdryVals(Surf, Mesh, Ty, 4, 'FY');

Wex(DofsFy) = Wex(DofsFy) + Fy;

% impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;

[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 3, 'UY');

BdryIdcs = DofsY;

disp([num2str(toc),'  Build matrices for OP']);

% define problem type and Nvar of optimization problem
if ( strcmp(stressState, 'PlaneStress') )       % Plane stress case
    NvarG = countGP;
    addNvar = 3*NvarG;
    Nvar = noDofs + NvarG + addNvar;
else
    NvarG = countGP;
    addNvar = 2*NvarG;
    Nvar = noDofs + NvarG + addNvar;
end

% define external work
Exf = zeros(Nvar,1);
for i = 1 : NvarG
    Exf(noDofs + i, 1) = sigmap*WeJ(i);
end

% boundary condition
sbc = length(BdryIdcs);
Bc = zeros(sbc, noDofs);
for i = 1 : sbc
    loca = BdryIdcs(i);
    Bc(i, loca) = 1;
end

Aeq = [Wex'; Bc];
beq = zeros(1 + sbc, 1);
beq(1,1) = 1; % external work of unitary

AeqIn = [];
% enforce incompressibility conditions
if ( strcmp(stressState, 'PlaneStrain') ) % Plane Strain case
    [AeqIn] = incompresibility(noDofs, Bstrain, countGP);
end
Aeq = [Aeq; AeqIn];
a1 = size(Aeq, 1);

% convert matrix constrains
A1 = [Aeq sparse(a1, NvarG + addNvar)];

% these variables appear when assigning C'k = r_i
A2 = addedVariables2D(noDofs, NvarG, addNvar, Bstrain, countGP, stressState);
aeq = [A1; A2];
b = [beq; zeros(addNvar + size(AeqIn, 1), 1)];

% cones for plane stress
if ( strcmp(stressState,'PlaneStress') )       % Plane Strain case
    for i = 1 : NvarG
        co(i, :) = [noDofs + i, noDofs+NvarG+3*i-2, noDofs+NvarG+3*i-1, noDofs+NvarG+3*i];
    end
else
    % cones for plane strain
    for i = 1 : NvarG
        co(i, :) = [noDofs + i, noDofs+NvarG+2*i-1, noDofs+NvarG+2*i];
    end
end

% call mosek optimisation tool
disp([num2str(toc),'  solving']);
clear prob
clear param
prob.c = Exf';
prob.a = aeq;
prob.buc = b';

prob.blc = b';
prob.blx = [];
prob.bux = [];

% define cones
prob.cell = cell(NvarG,1);
for i = 1:NvarG
    prob.cones{i}.type = 'MSK_CT_QUAD';
    prob.cones{i}.sub  = co(i,:);
end

param =[];
param.MSK_IPAR_PRESOLVE_USE = 'MSK_OFF';
[r,res] = mosekopt('minimize',prob,param);
try
    % display the optimal solution.
    u = res.sol.itr.xx;
    result = Exf'*u;
catch
    fprintf ('MSKERROR: Could not get solution')
end

loadfactor = result/sigmap

% export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 201), linspace(0, 1, 201)};
filename = 'SPGroovedRectPlateCollapse';
fieldname = 'Mechanisum';

d = u(1 : Mesh.NDof);
SPToVTK(Surf, d, ParaPts, filename, fieldname)

