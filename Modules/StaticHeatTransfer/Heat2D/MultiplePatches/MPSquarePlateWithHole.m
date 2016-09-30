% -----------------------------------------------------------
% Heat conduction in a square plate with a circular hole in the center
% BCs: prescribed temperature and heat flux
% this problem is modeled by four patches
% --------------------------------------------------------

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

format long
% control points
tic;
a = 0.3;
b = 1;
Surf = cell(1, 4);

% patch 1
CtrlPts{1} = zeros(4, 3, 2);

CtrlPts{1}(1 : 3, 1, 1) = [-b; -b; 0];
CtrlPts{1}(1 : 3, 2, 1) = [0; -b; 0];
CtrlPts{1}(1 : 3, 3, 1) = [b; -b; 0];

CtrlPts{1}(1 : 3, 1, 2) = [-a * cos(pi/4); -a * sin(pi/4); 0];
CtrlPts{1}(1 : 3, 2, 2) = [0; -a / cos(pi/4); 0];
CtrlPts{1}(1 : 3, 3, 2) = [a * cos(pi/4); -a * sin(pi/4); 0];

CtrlPts{1}(4, :, :) = 1;
CtrlPts{1}(:, 2, 2) = CtrlPts{1}(:, 2, 2) * cos(pi/4);

% knot vector

KntVect1{1} = [0 0 0 1 1 1];
KntVect1{2} = [0 0 1 1];

Surf{1} = CreateNURBS(KntVect1, CtrlPts{1});

% patch 2
CtrlPts{2} = zeros(4, 3, 2);

CtrlPts{2}(1 : 3, 1, 1) = [b; -b; 0];
CtrlPts{2}(1 : 3, 2, 1) = [b; 0; 0];
CtrlPts{2}(1 : 3, 3, 1) = [b; b; 0];

CtrlPts{2}(1 : 3, 1, 2) = [a * cos(pi/4); -a * sin(pi/4); 0];
CtrlPts{2}(1 : 3, 2, 2) = [a / cos(pi/4); 0; 0];
CtrlPts{2}(1 : 3, 3, 2) = [a * cos(pi/4); a * sin(pi/4); 0];

CtrlPts{2}(4, :, :) = 1;
CtrlPts{2}(:, 2, 2) = CtrlPts{2}(:, 2, 2) * cos(pi/4);
% knot vector

KntVect2{1} = [0 0 0 1 1 1];
KntVect2{2} = [0 0 1 1];

Surf{2} = CreateNURBS(KntVect2, CtrlPts{2});

% patch 3
CtrlPts{3} = zeros(4, 3, 2);

CtrlPts{3}(1 : 3, 1, 1) = [b; b; 0];
CtrlPts{3}(1 : 3, 2, 1) = [0; b; 0];
CtrlPts{3}(1 : 3, 3, 1) = [-b; b; 0];

CtrlPts{3}(1 : 3, 1, 2) = [a * cos(pi/4); a * sin(pi/4); 0];
CtrlPts{3}(1 : 3, 2, 2) = [0; a / cos(pi/4); 0];
CtrlPts{3}(1 : 3, 3, 2) = [-a * cos(pi/4); a * sin(pi/4); 0];

CtrlPts{3}(4, :, :) = 1;
CtrlPts{3}(:, 2, 2) = CtrlPts{3}(:, 2, 2) * cos(pi/4);
% knot vector

KntVect3{1} = [0 0 0 1 1 1];
KntVect3{2} = [0 0 1 1];

Surf{3} = CreateNURBS(KntVect3, CtrlPts{3});

% patch 4
CtrlPts{4} = zeros(4, 3, 2);

CtrlPts{4}(1 : 3, 1, 1) = [-b; b; 0];
CtrlPts{4}(1 : 3, 2, 1) = [-b; 0; 0];
CtrlPts{4}(1 : 3, 3, 1) = [-b; -b; 0];

CtrlPts{4}(1 : 3, 1, 2) = [-a * cos(pi/4); a * sin(pi/4); 0];
CtrlPts{4}(1 : 3, 2, 2) = [-a / cos(pi/4); 0; 0];
CtrlPts{4}(1 : 3, 3, 2) = [-a * cos(pi/4); -a * sin(pi/4); 0];

CtrlPts{4}(4, :, :) = 1;
CtrlPts{4}(:, 2, 2) = CtrlPts{4}(:, 2, 2) * cos(pi/4);
% knot vector

KntVect4{1} = [0 0 0 1 1 1];
KntVect4{2} = [0 0 1 1];

Surf{4} = CreateNURBS(KntVect4, CtrlPts{4});

% --------------------------------------------------------------

Interfaces = struct;
for i = 1 : 4
    Interfaces(i).Side1 = 2;
    Interfaces(i).Side2 = 1;
    Interfaces(i).Ornt = 1;
end
Interfaces(1).Patch1 = 1;
Interfaces(1).Patch2 = 2;

Interfaces(2).Patch1 = 2;
Interfaces(2).Patch2 = 3;

Interfaces(3).Patch1 = 3;
Interfaces(3).Patch2 = 4;

Interfaces(4).Patch1 = 4;
Interfaces(4).Patch2 = 1;

Boundaries = struct;

for i = 1 : 4
    Boundaries(i).Patches = i;
    Boundaries(i).Sides = 3;
end
Boundaries(5).Patches = [1; 2; 3; 4];
Boundaries(5).Sides = [4; 4; 4; 4];

Mesh = cell(1, numel(Surf));
for iPtc = 1 : numel(Surf)    
    Surf{iPtc} = KRefine(Surf{iPtc}, [6, 4], [4, 4], [0, 0]);
    Mesh{iPtc} = Mesh2D(Surf{iPtc}, 'ScalarField');
end

figure
hold on
grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
axis off
for i = 1 : numel(Surf)
    PlotGeo(Surf{i});
    PlotKnts(Surf{i});
%     PlotCtrlPts(Surf{i});
%     PlotCtrlNet(Surf{i});
end

[GNum, GDof] = MPInterfaceScalar2D(Interfaces, Mesh);

% material properties
k = 1; % W/(m * ^\circ C)
s = @(x, y) (2 * a / sqrt(x ^ 2 + y ^ 2) - 4);

NPatch = numel(Mesh);

% Assemble the system
disp([num2str(toc),'  Assembling the system'])

f = zeros(GDof, 1);
ii = cell(NPatch, 1);
jj = cell(NPatch, 1);
vv = cell(NPatch, 1);
for iPtc = 1 : NPatch
    [KVals, FVals] = calcLocalConductionMatrices2D(Mesh{iPtc}, Surf{iPtc}, k, s);
    [rsK, csK, vsK, FL] = convertToTripletStorage(Mesh{iPtc}, KVals, FVals);
    ii{iPtc} = GNum{iPtc}(rsK);
    jj{iPtc} = GNum{iPtc}(csK);
    vv{iPtc} = vsK;
    f(GNum{iPtc}) = f(GNum{iPtc}) + FL;
end
Rows = cell2mat(ii);
Cols = cell2mat(jj);
Vals = cell2mat(vv);

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
q0 = @(x, y) abs(2 * y * (1 - a / sqrt(x ^ 2 + y ^ 2)));

[flux, fluxDofs] = applyNewmannBdryVals(Surf, Mesh, q0, [1, 3], 'HFLUX', GNum, Boundaries);

f(fluxDofs) = f(fluxDofs) + flux;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
u0A = @(x, y) 0;
u0B = @(x, y) a ^ 2 + b ^ 2 + y ^ 2 - 2 * a * sqrt(y ^ 2 + b ^ 2);
[uA, BCDofsA] = projDrchltBdryVals(Surf, Mesh, u0A, 5, 'TEMP', GNum, Boundaries);
[uB, BCDofsB] = projDrchltBdryVals(Surf, Mesh, u0B, [2, 4], 'TEMP', GNum, Boundaries);

BdryIdcs = [BCDofsA; BCDofsB];
BdryVals = [uA; uB];

FreeIdcs = setdiff(1 : GDof, BdryIdcs);

d = zeros(GDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

% Export result to *.vtk file
disp([num2str(toc),'  Exporting result to *.vtk file'])
ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101)};

filename = 'MPSquarePlateWithHole';
fieldname = 'Temperature';
MPToVTK(Surf, d, GNum, ParaPts, filename, fieldname)

