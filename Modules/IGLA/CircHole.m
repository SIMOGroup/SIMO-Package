% Sheet with Small Circular Hole
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

close all
clear
clc

% define problem 
tic
a = 1;
L = 5;

% control points
CtrlPts1 = zeros(4, 4);
CtrlPts1(1 : 3, 1) = [0; L; 0];
CtrlPts1(1 : 3, 2) = [-L; L; 0];
CtrlPts1(1 : 3, 3) = [-L; L; 0];
CtrlPts1(1 : 3, 4) = [-L; 0; 0];

CtrlPts1(4, :) = 1;

Curv2 = NURBSCirc((a + L) / 2, [0, 0, 0], pi / 2, pi);
Curv2 = HRefine(Curv2, 1, 1/2);
CtrlPts2 = Curv2.CtrlPts4D;
clear Curv2

Curv3 = NURBSCirc(a, [0, 0, 0], pi / 2, pi);
Curv3 = HRefine(Curv3, 1, 1/2);
CtrlPts3 = Curv3.CtrlPts4D;
clear Curv3

KntVect{1} = [0 0 0 1/2 1 1 1];
KntVect{2} = [0 0 0 1 1 1];

CtrlPtsSurf = cat(3, CtrlPts1, CtrlPts2, CtrlPts3);
Surf = CreateNURBS(KntVect, CtrlPtsSurf);

% degree of basis function
p=2;q=2;
% repeated knot value inside knot vector
kx=2;ky=2;
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
% PlotCtrlPts(Surf);
PlotCtrlPts(Surf,1); % display control points
PlotCtrlNet(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'VectorField');

% assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])

% material properties
sigmap = 1; % Yield stress MPa
P1  = 1;     % load
P2  = 0;     % load

%  stressState ='PlaneStress';
stressState ='PlaneStrain';

% define problem type and total number of DOFs of system
if ( strcmp(stressState, 'PlaneStress') )       % Plane stress case
    noDofs=Mesh.NDof;
else
    imode=3;
    noDofs=Mesh.NDof+2*Mesh.NEl*imode;
end

[Bstrain, WeJ, countGP] = calcLocalBMatrices2D(Mesh, Surf, stressState);

Wex = zeros(noDofs, 1);

% impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])

Tx = @(x, y) -P1 * (abs(x + L) <= 1e-6);
Ty = @(x, y) P2 * (abs(y - L) <= 1e-6);

[Fx, DofsFx] = applyNewmannBdryVals(Surf, Mesh, Tx, 3, 'FX');
[Fy, DofsFy] = applyNewmannBdryVals(Surf, Mesh, Ty, 3, 'FY');

Wex(DofsFx) = Wex(DofsFx) + Fx;
Wex(DofsFy) = Wex(DofsFy) + Fy;

% impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;

[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 2, 'UY');
[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, h, 1, 'UX');

BdryIdcs = [DofsY; DofsX];

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

% call Mosek optimisation tool
disp([num2str(toc),'  solving']);
clear prob
clear param
prob.c = Exf';
prob.a = aeq;
prob.buc = b';
% axis off

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
    % display the optimal solution
    u = res.sol.itr.xx;
    result = Exf'*u;
catch
    fprintf ('MSKERROR: Could not get solution')
end

loadfactor = result/sigmap

% export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {setdiff(linspace(0, 1, 201), 0.5), linspace(0, 1, 201)};
filename = 'SPSheetCircHoleCollapse';
fieldname = 'Mechanisum';

d=u(1:Mesh.NDof);
SPToVTK(Surf, d, ParaPts, filename, fieldname)
