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
curv{1} = NURBSCirc(1, [0, 0, 0], 0, pi / 2);
curv{2} = NURBSCirc(2, [0, 0, 0], 0, pi / 2);

Surf = NURBSRuled(curv{Ri}, curv{Ro});
nel = 10;
p=2;q=4;
k1=1;k2=1;
Surf = KRefine(Surf, [nel, nel], [p, q], [p-k1, q-k2]);

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
E = 2e5;
nu = 0.3;

KVals = calcLocalStiffnessMatrices2D(Mesh, Surf, E, nu, 'PlaneStrain');
[Rows, Cols, Vals] = convertToTripletStorage(Mesh, KVals);

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

f = zeros(Mesh.NDof, 1);
% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
InnerPress = 6;
OuterPress = 0;

[Pi, DofsPi] = applyNewmannBdryVals(Surf, Mesh, @(x, y) InnerPress, 3, 'PRES');

f(DofsPi) = f(DofsPi) + Pi;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
h = @(x, y) 0;

[UY, DofsY] = projDrchltBdryVals(Surf, Mesh, h, 1, 'UY');
[UX, DofsX] = projDrchltBdryVals(Surf, Mesh, h, 2, 'UX');

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
% --------------------------------------------------------------------
% Evaluate IGA solution
[C, Stress] = StressEval(Surf, d, {linspace(0, 1, 401), 0.5}, E, nu,...
'PlaneStrain');

% Evaluate exact solution
x = reshape(C(1, :, :), 1, [])';
y = reshape(C(2, :, :), 1, [])';
alfa = atan2(y, x);
S = SPThickWalledCylindExactStress(Ri, Ro, InnerPress, OuterPress, nu);

figure
hold on
plot(rad2deg(alfa), S.xy(x, y), 'r', 'LineWidth', 1.5);
plot(rad2deg(alfa), reshape(squeeze(Stress(3, :, :)), [], 1), 'b')
legend('exact', 'approx')
% -----------------------------------------------------------------------

% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101)};
    
SPToVTKStress(Surf, d, ParaPts, E, nu, 'PlaneStrain',...
    'SPThickWalledCylindStress')

filename = 'SPThickWalledCylindDispl';
fieldname = 'Displ';
SPToVTK(Surf, d, ParaPts, filename, fieldname)

