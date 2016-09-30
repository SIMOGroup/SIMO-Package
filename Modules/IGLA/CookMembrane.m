% Trapezoidal panel 
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

% define problem
tic
L1=4.8;
L2=4.4;
L3=1.6;

% control points
CtrlPts = zeros(4, 2, 2);

CtrlPts(1 : 3, 1, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2, 1) = [L1; L2; 0];

CtrlPts(1 : 3, 1, 2) = [0; L2; 0];
CtrlPts(1 : 3, 2, 2) = [L1; L2+L3; 0];

CtrlPts(4, :, :) = 1;

KntVect{1} = [0 0 1 1];
KntVect{2} = [0 0 1 1];

% create NURBS surface in CAD 
Surf = CreateNURBS(KntVect, CtrlPts);

% degree of basis function
p=3;q=3;
% repeated knot value inside knot vector
kx=1;ky=1;
% number of elements in each direction
nelx=8; nely=8;
% h,p,k-refinements
Surf = KRefine(Surf, [nelx, nely], [p, q], [p-kx, q-ky]);

% show geometry
figure
hold on
axis equal
daspect([1, 1, 1])
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf,1);
%PlotCtrlPts(Surf);
PlotCtrlNet(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'VectorField');

% assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])
% material properties
sigmap = sqrt(3); % Yield stress MPa
P  = 1;     % load

stressState ='PlaneStress';
%  stressState ='PlaneStrain';

% define problem type and total number of DOFs of system
if ( strcmp(stressState, 'PlaneStress') )       % Plane stress case
    noDofs=Mesh.NDof;
else
    imode=3;
    noDofs=Mesh.NDof+2*Mesh.NEl*imode;
end

% tic
% [Bstrain, WeJ, countGP] = calcLocalBMatrices2D(Mesh, Surf, stressState);
% toc
% tic
[Bstrain, WeJ, countGP] = calcLocalBMatrices2DOp(Mesh, Surf, stressState);
% toc

Wex = zeros(noDofs, 1);

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])

Ty = @(x, y) -P;

[Fy1, DofsFy1] = applyNewmannBdryVals(Surf, Mesh, Ty, 2, 'FY');

Wex(DofsFy1) = Wex(DofsFy1) + Fy1;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;

[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, h, 1, 'UX');
[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 1, 'UY');

BdryIdcs = [DofsX; DofsY];

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

loadfactor = result

% export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 100), linspace(0, 1, 100)};

%     
% SPToVTKStress(Surf, d, ParaPts, E, nu, 'PlaneStress',...
%     'Trapezoidalpanel')
% 
 filename = 'TrapezoidalpanelCollapseOp';
 fieldname = 'Mechanisum';
 d=u(1:Mesh.NDof);
 SPToVTK(Surf, d, ParaPts, filename, fieldname)
% 
