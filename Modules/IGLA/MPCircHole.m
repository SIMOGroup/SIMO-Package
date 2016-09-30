% Sheet with Small Circular Hole using Multi patches

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
R = 1;
L = 5;

% patch 1
CtrlPts{1} = zeros(4, 2, 3);

CtrlPts{1}(1 : 3, 1, 1) = [-L; 0; 0];
CtrlPts{1}(1 : 3, 2, 1) = [-R; 0; 0];

CtrlPts{1}(1 : 3, 1, 2) = [-L; L / 2; 0];
CtrlPts{1}(1 : 3, 2, 2) = [-R; R * tan(pi / 8); 0];

CtrlPts{1}(1 : 3, 1, 3) = [-L; L; 0];
CtrlPts{1}(1 : 3, 2, 3) = [-R * cos(pi / 4); R * sin(pi / 4); 0];

CtrlPts{1}(4, :, :) = 1;

% patch 2
CtrlPts{2} = zeros(4, 2, 3);

CtrlPts{2}(1 : 3, 1, 1) = [-L; L; 0];
CtrlPts{2}(1 : 3, 2, 1) = [-R * cos(pi / 4); R * sin(pi / 4); 0];

CtrlPts{2}(1 : 3, 1, 2) = [-L / 2; L; 0];
CtrlPts{2}(1 : 3, 2, 2) = [-R * tan(pi / 8); R; 0];

CtrlPts{2}(1 : 3, 1, 3) = [0; L; 0];
CtrlPts{2}(1 : 3, 2, 3) = [0; R; 0];

CtrlPts{2}(4, :, :) = 1;

W = cos(pi / 8);

CtrlPts{1}(:, 2, 2) = CtrlPts{1}(:, 2, 2) * W;
CtrlPts{2}(:, 2, 2) = CtrlPts{2}(:, 2, 2) * W;

KntVect{1} = [0 0 1 1];
KntVect{2} = [0 0 0 1 1 1];

Surf{1} = CreateNURBS(KntVect, CtrlPts{1});
Surf{2} = CreateNURBS(KntVect, CtrlPts{2});

% ---------------------------------------------------

Mesh = cell(1, numel(Surf));
for iPtc = 1 : numel(Surf)    
    Surf{iPtc} = KRefine(Surf{iPtc}, [5, 5], [2, 2], [0, 0]);
    Mesh{iPtc} = Mesh2D(Surf{iPtc}, 'VectorField');
end

figure
hold on
axis equal
daspect([1, 1, 1])

for i = 1 : numel(Surf)
    PlotGeo(Surf{i});
    PlotKnts(Surf{i});
    PlotCtrlPts(Surf{i});
    PlotCtrlNet(Surf{i});
end
% --------------------------------------------------------------
Interfaces = struct;
Interfaces(1).Patch1 = 1;
Interfaces(1).Patch2 = 2;
Interfaces(1).Side1 = 4;
Interfaces(1).Side2 = 3;
Interfaces(1).Ornt = 1;

Boundaries = struct;

Boundaries(1).NSides = 1;
Boundaries(1).Patches = 1;
Boundaries(1).Sides = 1;

Boundaries(2).NSides = 1;
Boundaries(2).Patches = 1;
Boundaries(2).Sides = 3;

Boundaries(3).NSides = 1;
Boundaries(3).Patches = 2;
Boundaries(3).Sides = 4;

% Set a local-to-global numbering for all patches
[GNum, GDof] = MPInterfaceVector2D(Interfaces, Mesh);

NPatch = numel(Mesh);

% material properties
sigmap = 1; % Yield stress MPa
P1 = 1;
stressState ='PlaneStress';

[Bstrain, WeJ, countGP] = calcLocalBMatrices2D_mod(Mesh, Surf, GNum, GDof);

f = zeros(GDof,1);

% Impose natural boundary conditions
disp([num2str(toc),'Imposing natural boundary conditions'])

Tx = @(x, y) -P1;
[Fx1, DofsFx1] = applyNewmannBdryVals(Surf, Mesh, Tx, 1, 'FX', GNum, Boundaries);
f(DofsFx1) = f(DofsFx1) + Fx1;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;
[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 2, 'UY', GNum, Boundaries);
[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, h, 3, 'UX', GNum, Boundaries);

BdryIdcs = [DofsY; DofsX];

disp([num2str(toc),'  Build matrices for OP']);
% +----------------------------+
% | OBJECTIVE FUNCTION AMPHA+  |
% +----------------------------+
if ( strcmp(stressState, 'PlaneStress') )       % Plane stress case
    NvarG = countGP;
    addNvar = 3*NvarG;
else
    NvarG = countGP;
    addNvar = 2*NvarG;
end
Nvar = GDof + NvarG + addNvar;

% External loading
Wex = f';
clear f;

Exf = zeros(Nvar,1);
for i = 1 : NvarG
    Exf(GDof + i, 1) = sigmap*WeJ(i);
end

% boundary condition
sbc = length(BdryIdcs);
Bc = zeros(sbc, GDof);
for i = 1 : sbc
    loca = BdryIdcs(i);
    Bc(i, loca) = 1;
end

Aeq = [Wex; Bc];
beq = zeros(1 + sbc, 1);
beq(1, 1) = 1; % external work of unitary

AeqIn = [];
% Enforce incompressibility conditions
if ( strcmp(stressState, 'PlaneStrain') ) % Plane Strain case
    [AeqIn] = incompresibility(GDof, Bstrain, countGP);
end
Aeq = [Aeq; AeqIn];
a1 = size(Aeq, 1);

% trans to Xadded_variables_plate.m
A1 = [Aeq, sparse(a1, NvarG + addNvar)];

% these variables appear when assigning C'k = r_i
A2 = addedVariables2D(GDof, NvarG, addNvar, Bstrain, countGP, stressState);
aeq = [A1; A2];
b = [beq; zeros(addNvar + size(AeqIn, 1), 1)];

% cones for plane stress
if ( strcmp(stressState,'PlaneStress') )       % Plane Stress case
    for i = 1 : NvarG
        co(i, :) = [GDof + i, GDof + NvarG + 3*i - 2, GDof + NvarG + 3*i - 1, GDof + NvarG + 3*i];
    end
else
    % cones for plane strain
    for i = 1 : NvarG
        co(i, :) = [GDof + i, GDof + NvarG + 2*i-1, GDof + NvarG + 2*i];
    end
end

% +--------------+
% | OPTIMIZATION |
% +--------------+
% Mosek optimisation tool

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
    prob.cones{i}.sub  = co(i, :);
end

param =[];
param.MSK_IPAR_PRESOLVE_USE = 'MSK_OFF';
[r,res] = mosekopt('minimize',prob,param);
try
    % Display the optimal solution.
    u = res.sol.itr.xx;
    result = Exf'*u;
catch
    fprintf ('MSKERROR: Could not get solution')
end

loadfactor = result/sigmap

d = u(1 : GDof);
% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101)};

MPToVTK(Surf, d, GNum, ParaPts, 'MPCircHoleCollapse', 'Mechanisum')

