clear 
close all
clc

tic
curv{1} = NURBSCirc(1, [0, 0, 0], 0, pi / 2);
curv{2} = NURBSCirc(2, [0, 0, 0], 0, pi / 2);

Surf = NURBSRuled(curv{1}, curv{2});

Volu = NURBSExtrude(Surf, [0 0 2]);

Volu = KRefine(Volu, [8, 4, 10], [4 4 4], [3 3 3]);
% Volu = KRefine(Volu, [2, 2, 2], [2 2 2], [1 1 1]);

figure
hold on
grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

PlotGeo(Volu)
PlotKnts(Volu)
PlotCtrlPts(Volu)
PlotCtrlNet(Volu)

Mesh = Mesh3D(Volu, 'VectorField');

% Assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])

E = 2e11;  % Young modulus
nu = 0.3;  % Poissons ratio

KVals = calcLocalStiffnessMatrices3D(Mesh, Volu, E, nu);
[Rows, Cols, Vals] = convertToTripletStorage(Mesh, KVals);

% Convert triplet data to sparse matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

press = @(x, y, z) 1;

f = zeros(Mesh.NDof, 1);

% Impose natural boundary conditions
disp([num2str(toc),'  Imposing natural boundary conditions'])
[Press, PressDofs] = applyNewmannBdryVals(Volu, Mesh, press, 3, 'PRES');

f(PressDofs) = f(PressDofs) + Press;

symmDofs = unique([Mesh.Boundary(1).CompDofs{2},...
    Mesh.Boundary(5).CompDofs{3}, Mesh.Boundary(2).CompDofs{1},...
    Mesh.Boundary(6).CompDofs{3}]);

BdryIdcs = symmDofs;
BdryVals = zeros(1, numel(BdryIdcs))';

FreeIdcs = setdiff(1 : Mesh.NDof, BdryIdcs);

d = zeros(Mesh.NDof, 1);
d(BdryIdcs) = BdryVals;

f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

% Export solution to vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 21), linspace(0, 1, 11), linspace(0, 1, 11)};

filename = 'SPSolidCylindDispl';
fieldname = 'Displ';
SPToVTK(Volu, d, ParaPts, filename, fieldname)

filename1 = 'SPSolidCylindStress';
SPToVTKStress(Volu, d, ParaPts, E, nu, '3D', filename1);






