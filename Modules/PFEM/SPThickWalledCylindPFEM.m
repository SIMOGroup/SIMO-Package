% Thick-walled cylinder

close all
clear
clc
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points
tic;
Ri = 1;
Ro = 2;

CtrlPts = zeros(4, 2, 3);
CtrlPts(4, :, :) = 1;
CtrlPts(1 : 3, 1, 1) = [Ri; 0; 0];
CtrlPts(1 : 3, 2, 1) = [Ro; 0; 0];

CtrlPts(1 : 3, 1, 2) = [Ri; Ri; 0];
CtrlPts(1 : 3, 2, 2) = [Ro; Ro; 0];

CtrlPts(1 : 3, 1, 3) = [0; Ri; 0];
CtrlPts(1 : 3, 2, 3) = [0; Ro; 0];

fac = 1 / sqrt(2);
CtrlPts(:, :, 2) = CtrlPts(:, :, 2) * fac;
KntVect = {[0, 0, 1, 1], [0, 0, 0, 1, 1, 1]};
Surf = CreateNURBS(KntVect, CtrlPts);
Surf = KRefine(Surf, [8, 16], [2, 2], [1, 1]);

figure
hold on
axis equal
daspect([1, 1, 1])

PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf, 1);
% PlotCtrlNet(Surf);

% ---------------------------------------------------------------------
Mesh = Mesh2D(Surf, 'VectorField');

% Assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])
% material properties
E = 21000;
nu = 0.3;
tic
KVals = calcLocalStiffnessMatrices2D(Mesh, Surf, E, nu, 'PlaneStrain');
toc
% tic
% KVals = calcLocalStiffnessMatrices2DPFEM(Mesh, Surf, E, nu, 'PlaneStrain');
% toc

[Rows, Cols, Vals] = convertToTripletStorage(Mesh, KVals);
clear KVals

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);

clear Rows Cols Vals 

% spy(K)
f = zeros(Mesh.NDof, 1);
% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
InnerPress = -6;
OuterPress = 0;

[Vals, GDofs] = applyNewmannBdryVals(Surf, Mesh,  @(x, y) InnerPress, 1, 'PRES');

f(GDofs) = f(GDofs) + Vals;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;

[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 3, 'UY');
[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, h, 4, 'UX');

BdryIdcs = [DofsY; DofsX];
BdryVals = [UY; UX];

FreeIdcs = setdiff(1 : Mesh.NDof, BdryIdcs);

d = zeros(Mesh.NDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);


% Compare stress between IGA and exact solution
% ---------------------------------------------------------------------
% Evaluate IGA solution
[C, Stress] = StressEval(Surf, d, {0.5, linspace(0, 1, 401)}, E, nu,...
'PlaneStrain');

% Evaluate exact solution
x = reshape(C(1, :, :), 1, [])';
y = reshape(C(2, :, :), 1, [])';
alfa = atan2(y, x);
S = SPThickWalledCylindExactStress(Ri, Ro, -InnerPress, OuterPress, nu);

figure
hold on
plot(rad2deg(alfa), S.xy(x, y), 'r', 'LineWidth', 1.5);
plot(rad2deg(alfa), reshape(squeeze(Stress(3, :, :)), [], 1), 'b')
legend('exact', 'approx')
% ---------------------------------------------------------------------

StrainEnergy = 0.5 * f' * d * 4

% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101)};
    
SPToVTKStress(Surf, d, ParaPts, E, nu, 'PlaneStrain', 'SPThickWalledCylindStressSmooth')
SPToVTK(Surf, d, ParaPts, 'SPThickWalledCylindDisplSmooth', 'Displ')

