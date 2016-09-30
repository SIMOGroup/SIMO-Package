% ------------------------------------------------------------------------
clear 
close all
clc 
addpath('C:\Program Files\Mosek\7\toolbox\r2013aom');

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

tic
R = 1; % m
L = 2 * R;
%----------------- Patch 01 --------------
Curv1  = NURBSCirc(R, [0, 0, 0], 0, pi / 4);
Curv2  = NURBSLine([L, 0], [L, L]);
Surf{1} = NURBSRuled(Curv1, Curv2);
%----------------- Patch 02 --------------
Curv1  = NURBSCirc(R, [0, 0, 0], pi / 4, pi / 2);
Curv2  = NURBSLine([L, L], [0, L]);
Surf{2} = NURBSRuled(Curv1, Curv2);
%----------------- Patch 03 --------------
Curv1  = NURBSCirc(R, [2*L, 0], 3*pi / 4, pi);
Curv2  = NURBSLine([L, L], [L, 0]);
Surf{3} = NURBSRuled(Curv1, Curv2);
%----------------- Patch 04 --------------
Curv1  = NURBSCirc(R, [2*L, 0], pi / 2, 3*pi / 4);
Curv2  = NURBSLine([2*L, L], [L, L]); 
Surf{4} = NURBSRuled(Curv1, Curv2);
%----------------- Patch 05 --------------
Curv1  = NURBSLine([L, L], [0, L]);
Curv2  = NURBSLine([L, 2*L],[0, 2*L]);
Surf{5} = NURBSRuled(Curv1, Curv2);
%----------------- Patch 06 --------------
Curv1  = NURBSLine([2*L, L], [L, L]);
Curv2  = NURBSLine(2*[L, L], [L, 2*L]);
Surf{6} = NURBSRuled(Curv1, Curv2);

Interfaces = struct;
Side1  =[2, 4, 1, 4, 4, 1];
Side2  =[1, 4, 2, 3, 3, 2];
Patch1 =[1, 1, 3, 2, 4, 5];
Patch2 =[2, 3, 4, 5, 6, 6];
Ornt   =[1, -1, 1, 1, 1, 1];
for iConect = 1 : numel(Side1)
    Interfaces(iConect).Side1 = Side1(iConect);
    Interfaces(iConect).Side2 = Side2(iConect);
    Interfaces(iConect).Ornt = Ornt(iConect);
    Interfaces(iConect).Patch1 = Patch1(iConect);
    Interfaces(iConect).Patch2 = Patch2(iConect);
end
%--------------------------------------------------------------
Boundaries = struct;
Boundaries(1).Patches = [1; 3];
Boundaries(1).Sides = [1; 2];
Boundaries(2).Patches = [3; 4];
Boundaries(2).Sides = [3; 3];
Boundaries(3).Patches = [4; 6];
Boundaries(3).Sides = [1; 1];
Boundaries(4).Patches = [6; 5];
Boundaries(4).Sides = [4; 4];
Boundaries(5).Patches = [5; 2];
Boundaries(5).Sides = [2; 2];
Boundaries(6).Patches = [2; 1];
Boundaries(6).Sides = [3; 3];
%---------------------------------------------------------------

Mesh = cell(numel(Surf), 1);
for iPtc = 1:numel(Surf)
    Surf{iPtc} = KRefine(Surf{iPtc}, [5, 5], [2, 2], [0, 0]);
    Mesh{iPtc} = Mesh2D(Surf{iPtc}, 'VectorField');
end
% Set a local-to-global numbering for all patches
[GNum, GDof] = MPInterfaceVector2D(Interfaces, Mesh);

NPatch = numel(Mesh);

figure
hold on
axis equal
daspect([1, 1, 1])
for iPtc = 1:numel(Surf)
    PlotGeo(Surf{iPtc});
    PlotKnts(Surf{iPtc});
    PlotCtrlPts(Surf{iPtc}, 1, GNum{iPtc});
end
% material properties
sigmap = 1; % Yield stress MPa
P1 = 1;
stressState ='PlaneStress';

%[Bstrain, WeJ, countGP] = calcLocalBMatrices2D_mod(Mesh, Surf, GNum, GDof);

[Bstrain, WeJ, countGP] = calcLocalBMatrices2DOp(Mesh, Surf, stressState);

f = zeros(GDof,1);

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])

Ty = @(x, y) P1;
[Fy1, DofsFy1] = applyNewmannBdryVals(Surf, Mesh, Ty, 4, 'FY', GNum, Boundaries);
f(DofsFy1) = f(DofsFy1) + Fy1;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;
[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 1, 'UY', GNum, Boundaries);

BdryIdcs = DofsY;

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
    prob.cones{i}.sub  = co(i,:);
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

MPToVTK(Surf, d, GNum, ParaPts, 'MPGroovedRectPlateCollapse', 'Mechanisum')
