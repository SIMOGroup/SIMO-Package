% L-shape 3 patches
clear
close all
clc

% control points
tic;
L = 1;
D = 0.5;

CtrlPts = cell(3, 1);
% patch 1
CtrlPts{1} = zeros(4, 2, 2);

CtrlPts{1}(1 : 3, 1, 1) = [0; L; 0];
CtrlPts{1}(1 : 3, 2, 1) = [0; D; 0];

CtrlPts{1}(1 : 3, 1, 2) = [D; L; 0];
CtrlPts{1}(1 : 3, 2, 2) = [D; D; 0];

CtrlPts{1}(4, :, :) = 1;

% knot vector

KntVect1 = [0 0 1 1];
KntVect2 = [0 0 1 1];

Surf{1} = CreateNURBS({KntVect1, KntVect2}, CtrlPts{1});

% patch 2
CtrlPts{2} = zeros(4, 2, 2);

CtrlPts{2}(1 : 3, 1, 1) = [0; D; 0];
CtrlPts{2}(1 : 3, 2, 1) = [0; 0; 0];

CtrlPts{2}(1 : 3, 1, 2) = [D; D; 0];
CtrlPts{2}(1 : 3, 2, 2) = [D; 0; 0];

CtrlPts{2}(4, :, :) = 1;

% knot vector

KntVect1 = [0 0 1 1];
KntVect2 = [0 0 1 1];

Surf{2} = CreateNURBS({KntVect1, KntVect2}, CtrlPts{2});

% patch 3
CtrlPts{3} = zeros(4, 2, 2);

CtrlPts{3}(1 : 3, 1, 1) = [D; 0; 0];
CtrlPts{3}(1 : 3, 2, 1) = [L; 0; 0];

CtrlPts{3}(1 : 3, 1, 2) = [D; D; 0];
CtrlPts{3}(1 : 3, 2, 2) = [L; D; 0];

CtrlPts{3}(4, :, :) = 1;

% knot vector

KntVect1 = [0 0 1 1];
KntVect2 = [0 0 1 1];

Surf{3} = CreateNURBS({KntVect1, KntVect2}, CtrlPts{3});

figure
hold on
axis equal
 for i = 1 : numel(Surf)
     daspect([1, 1, 1])
     PlotGeo(Surf{i});
     PlotKnts(Surf{i});
     PlotCtrlPts(Surf{i},1);
     PlotCtrlNet(Surf{i});
 end

Interfaces = struct;
Side1  = [2, 4];
Side2  = [1, 1];
Patch1 = [1, 2];
Patch2 = [2, 3];

Ornt =[1, 1];
for iConect = 1 : numel(Side1)
    Interfaces(iConect).Side1 = Side1(iConect);
    Interfaces(iConect).Side2 = Side2(iConect);
    Interfaces(iConect).Ornt = Ornt(iConect);
    Interfaces(iConect).Patch1 = Patch1(iConect);
    Interfaces(iConect).Patch2 = Patch2(iConect);
end

%--------------------------------------------------------------

Boundaries = struct;
Boundaries(1).Patches = [1; 2];
Boundaries(1).Sides = [3; 3];
Boundaries(2).Patches = [2; 3];
Boundaries(2).Sides = [2; 3];
Boundaries(3).Patches = 3;
Boundaries(3).Sides = 2;
Boundaries(4).Patches = 3;
Boundaries(4).Sides = 4;
Boundaries(5).Patches = 1;
Boundaries(5).Sides = 4;
Boundaries(6).Patches = 1;
Boundaries(6).Sides = 1;

Mesh = cell(numel(Surf), 1);

for iPtc = 1:numel(Surf)
%Surf{iPtc} = KRefine(Surf{iPtc}, [5, 5], [2, 2], [0, 0]);
    Mesh{iPtc} = Mesh2D(Surf{iPtc}, 'VectorField');
end

%Set a local-to-global numbering for all patches
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

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
%q0 = @(x, y) abs(2 * y * (1 - a / sqrt(x ^ 2 + y ^ 2)));
q0 = 

[flux, fluxDofs] =  applyNewmannBdryVals(Surf, Mesh, q0, [1, 3], 'HFLUX', GNum, Boundaries);

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


