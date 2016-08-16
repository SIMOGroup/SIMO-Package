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

CtrlPts{1}(1 : 3, 1, 1) = [0; 0; 0];
CtrlPts{1}(1 : 3, 2, 1) = [D; 0; 0];

CtrlPts{1}(1 : 3, 1, 2) = [0; D; 0];
CtrlPts{1}(1 : 3, 2, 2) = [D; D; 0];

CtrlPts{1}(4, :, :) = 1;

% knot vector

KntVect1 = [0 0 1 1];
KntVect2 = [0 0 1 1];

Surf{1} = CreateNURBS({KntVect1, KntVect2}, CtrlPts{1});

% patch 2
CtrlPts{2} = zeros(4, 2, 2);

CtrlPts{2}(1 : 3, 1, 1) = [0; D; 0];
CtrlPts{2}(1 : 3, 2, 1) = [D; D; 0];

CtrlPts{2}(1 : 3, 1, 2) = [0; L; 0];
CtrlPts{2}(1 : 3, 2, 2) = [D; L; 0];

CtrlPts{2}(4, :, :) = 1;

% knot vector

KntVect1 = [0 0 1 1];
KntVect2 = [0 0 1 1];

Surf{2} = CreateNURBS({KntVect1, KntVect2}, CtrlPts{2});

% patch 3
CtrlPts{3} = zeros(4, 2, 2);

CtrlPts{3}(1 : 3, 1, 1) = [D; D; 0];
CtrlPts{3}(1 : 3, 2, 1) = [L; D; 0];

CtrlPts{3}(1 : 3, 1, 2) = [D; L; 0];
CtrlPts{3}(1 : 3, 2, 2) = [L; L; 0];

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
Side1  = [4, 2];
Side2  = [3, 1];
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
Boundaries(1).Sides = [1; 1];
Boundaries(2).Patches = 1;
Boundaries(2).Sides = 3;
Boundaries(3).Patches = 1;
Boundaries(3).Sides = 2;
Boundaries(4).Patches = 3;
Boundaries(4).Sides = 3;
Boundaries(5).Patches = 3;
Boundaries(5).Sides = 2;
Boundaries(6).Patches = [3; 2];
Boundaries(6).Sides = [4,4];

Mesh = cell(numel(Surf), 1);

for iPtc = 1:numel(Surf)
    Surf{iPtc} = KRefine(Surf{iPtc}, [5, 5], [2, 2], [1, 1]);
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
% Material properties
E = 10000;
nu = 0.3;
internalLength = 0.01;

% Assemble the system
disp([num2str(toc),'  Assembling the system'])

f = zeros(GDof, 1);
ii = cell(NPatch, 1);
jj = cell(NPatch, 1);
vv = cell(NPatch, 1);
for iPtc = 1 : NPatch
    %[KVals, FVals] = calcLocalStiffnessMatrices2D_Gradient(Mesh{iPtc}, Surf{iPtc}, k, s);
    KVals = calcLocalStiffnessMatrices2D_Gradient(Mesh{iPtc}, Surf{iPtc}, E, nu, internalLength, 'PlaneStrain');
    %[rsK, csK, vsK, FL] = convertToTripletStorage(Mesh{iPtc}, KVals, FVals);
    [rsK, csK, vsK] = convertToTripletStorage(Mesh{iPtc}, KVals);
    ii{iPtc} = GNum{iPtc}(rsK);
    jj{iPtc} = GNum{iPtc}(csK);
    vv{iPtc} = vsK;
    %f(GNum{iPtc}) = f(GNum{iPtc}) + FL;
end
Rows = cell2mat(ii);
Cols = cell2mat(jj);
Vals = cell2mat(vv);
% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
%q0 = @(x, y) abs(2 * y * (1 - a / sqrt(x ^ 2 + y ^ 2)));
q0 = @(x, y) -1;

%[flux, fluxDofs] = applyNewmannBdryVals(Surf, Mesh, q0, [1, 3], 'HFLUX', GNum, Boundaries);
[Fx, DofsFx] = applyNewmannBdryVals(Surf, Mesh, q0, 1, 'FX', GNum, Boundaries);

%f(fluxDofs) = f(fluxDofs) + flux;
f(DofsFx) = f(DofsFx) + Fx;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
%u0A = @(x, y) 0;
%u0B = @(x, y) a ^ 2 + b ^ 2 + y ^ 2 - 2 * a * sqrt(y ^ 2 + b ^ 2);
%[uA, BCDofsA] = projDrchltBdryVals(Surf, Mesh, u0A, 5, 'TEMP', GNum, Boundaries);
%[uB, BCDofsB] = projDrchltBdryVals(Surf, Mesh, u0B, [2, 4], 'TEMP', GNum, Boundaries);

%BdryIdcs = [BCDofsA; BCDofsB];
%BdryVals = [uA; uB];

h = @(x, y) 0;

[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 2, 'UY', GNum, Boundaries);
[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, h, 5, 'UX', GNum, Boundaries);

BdryIdcs = [DofsY; DofsX];
BdryVals = [UY; UX];

FreeIdcs = setdiff(1 : GDof, BdryIdcs);

d = zeros(GDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

% Evaluate IGA solution
%[C, Stress] = StressEval(Surf, d, {linspace(0, 1, 401), 0.5}, E, nu,'PlaneStrain');

% Export result to *.vtk file
disp([num2str(toc),'  Exporting result to *.vtk file'])
ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101)};

filename = 'Lplate-Multipatch';
fieldname = 'Diplacement';
MPToVTK(Surf, d, GNum, ParaPts, filename, fieldname)

% for iPtc = 1 : NPatch
    ParaPts = {linspace(0, 1, 201), linspace(0, 1, 201)};
    MPToVTKStress(Surf, d, GNum, ParaPts, E, nu, 'PlaneStrain', 'Lplate')
    %SPToVTK(Surf, d, ParaPts, 'CantiliverBeamDisplGrad', 'Displ')
% end





