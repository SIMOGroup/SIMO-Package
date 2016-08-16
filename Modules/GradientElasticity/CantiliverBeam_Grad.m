% CantiliverBeam - Gradient Elasticity 
% ref: [Harm Askes, Elias C. Aifantis]_Numerical modeling of size effects with gradient elasticity - Formulation, meshless discretization and examples

close all
clear
clc
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points
tic;
L = 36;
D = 12;
CtrlPts = zeros(4, 2, 2);

CtrlPts(1 : 3, 1, 1) = [0; -D/2; 0];
CtrlPts(1 : 3, 2, 1) = [L; -D/2; 0];

CtrlPts(1 : 3, 1, 2) = [0; D/2; 0];
CtrlPts(1 : 3, 2, 2) = [L; D/2; 0];  %Harm Askes

CtrlPts(4, :, :) = 1;

KntVect{1} = [0 0 1 1];
KntVect{2} = [0 0 1 1];

% degree of basis function
p = 2; q = 2;
% repeated knot value inside knot vector
kx = 1; ky = 1;
% number of elements in each direction
%nelx = 30; nely = 10;
nelx = 2; nely = 2;
% create NURBS surface in CAD 
Surf = CreateNURBS(KntVect, CtrlPts);
% h,p,k-refinements
Surf = KRefine(Surf, [nelx, nely], [p, q], [p-kx, q-ky]);

figure
hold on
axis equal
daspect([1, 1, 1])
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf,1);
PlotCtrlNet(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'VectorField');

% Assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])
% material properties
E = 1e5;
nu = 0.25;

internalLength = linspace(D, D / 200, 20);
wh = zeros(2, numel(internalLength));

%Iteration

for l = 1 : numel(internalLength)
    KVals = calcLocalStiffnessMatrices2D_Gradient(Mesh, Surf, E, nu, internalLength(l), 'PlaneStrain');

    [Rows, Cols, Vals] = convertToTripletStorage(Mesh, KVals);

    % Convert triplet data to sparse matrix
    K = sparse(Rows, Cols, Vals);
    clear Rows Cols Vals

    f = zeros(Mesh.NDof, 1);
    % Impose natural boundary conditions
    disp([num2str(toc),'  Imposing natural boundary conditions'])
    P = 1000;
    I = D ^ 3 / 12;

    Ty = @(x, y) -P * (D ^ 2 / 4 - y ^ 2) / (2 * I);

    [Fy, DofsFy] = applyNewmannBdryVals(Surf, Mesh, Ty, 2, 'FY');
    f(DofsFy) = f(DofsFy) + Fy;

    % Impose essential boundary conditions
    disp([num2str(toc),'  Imposing essential boundary conditions'])

    % constrain displacement
    [UX, DofsX] = projDrchltBdryVals(Surf, Mesh, @(x, y) 0, 1, 'UX');
    [UY, DofsY] = projDrchltBdryVals(Surf, Mesh, @(x, y) 0, 1, 'UY');

    % constrain rotation
    %adjacentDofs = Mesh.Boundary(1).NextLayerDofs.CompDofs{2};
    adjacentDofs =[]

    BdryIdcs = [DofsY; DofsX; adjacentDofs];
    BdryVals = [UY; UX; zeros(numel(adjacentDofs), 1)];

    [f, d, FreeIdcs] = applyDrchltBdryVals(BdryIdcs, BdryVals, K, f); 

    % Solve the system
    disp([num2str(toc),'  Solving the system'])
    d(FreeIdcs) = K(FreeIdcs, FreeIdcs) \ f(FreeIdcs);

    % deflection at tip

    [C, wh(:, l)] = NURBSEval(Surf, {1, 1}, d);
end

ka = E * D ^ 3 / (4 * L ^ 3);   % analytical stiffness
kh = P ./ abs(wh(2, :));        % numerical stiffness

figure
hold on
grid on
plot(log(D ./ internalLength), log(ka ./ kh))

% % Export solution to vtk format
% disp([num2str(toc),'  Exporting result to *.vtk format'])
% ParaPts = {linspace(0, 1, 201), linspace(0, 1, 201)};
% SPToVTKStress(Surf, d, ParaPts, E, nu, 'PlaneStress', 'CantiliverBeamStressGrad')
% SPToVTK(Surf, d, ParaPts, 'CantiliverBeamDisplGrad', 'Displ')

