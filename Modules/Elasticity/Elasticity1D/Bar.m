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

clc
clear
close all

tic;

% material properties
E = @(x) 1; % Young's Modulus

% physical parameters
A = @(x) 1; % Area of the bar
L = 1; % Length of the bar

% essential boundary condition
u0 = 0.5;

% natural boundary conditions
q1 = -3;
q2 = 5;
P = -1;
% --------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points

CtrlPts = zeros(4, 2);

CtrlPts(1 : 3, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2) = [L; 0; 0];

CtrlPts(4, :) = 1;

% knot vector

KntVect{1} = [0 0 1 1];
p = 3;
nel = 1;
k = 1;
Line = CreateNURBS(KntVect,CtrlPts);
Line = KRefine(Line, nel, p, p - k); % use for h,p,k-refinement

figure
hold on
axis equal
daspect([1, 1, 1])
PlotGeo(Line);
PlotKnts(Line);
PlotCtrlPts(Line);
PlotCtrlNet(Line);

Mesh = Mesh1D(Line);
% -------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------PROCESSING--------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'  Assembling the system'])
df = @(x) q1 * (1 - x / L) + q2 * x / L;

[KVals, FVals] = calcLocalStiffnessMatrices1D(Mesh, Line, E, A, df);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals);

% assemble elementary stiffness matrix to the global matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

% Impose Dirichlet boundary conditions
disp([num2str(toc),'  Imposing Dirichlet boundary conditions'])
BCIdx = 1;
BCVal = 0.5;
FreeIdx = setdiff(1 : Mesh.NDof, BCIdx);
d = zeros(Mesh.NDof, 1);
d(BCIdx) = BCVal;
f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdx) * BCVal; 

% Impose Neumann boundary conditions
disp([num2str(toc),'  Imposing Neumann boundary conditions'])
f(end) = f(end) + P;

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdx) = K(FreeIdx, FreeIdx) \ f(FreeIdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------- POST-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot control points
figure
hold on
title('Displacement')
plot(Line.CtrlPts3D(1, :), d, 'r.','MarkerSize', 20);
plot(Line.CtrlPts3D(1, :), d, 'k--');
set(gcf, 'color', 'white');
% interpolate solution field
ParaPts = linspace(0, 1, 1001);
[C, U] = NURBSEval(Line, {ParaPts}, d);    
ApproxDispl = plot(C(1, :), U, 'k-', 'LineWidth', 1);
x = C(1, :);

% Evaluate exact solution
solu = BarExactSolution(q1, q2, E, A, L, u0, P);
ExactDispl = plot(x, solu.displ(x),'r-', 'LineWidth', 1);
legend([ApproxDispl, ExactDispl], 'approx', 'exact')

figure
hold on
title('Stress')
set(gcf, 'color', 'white');
g = GradEval(Line.KntVect, Line.CtrlPts4D, {ParaPts}, d'); % evaluate strain
ApproxStress = plot(C(1, :), g * E(x), 'k--','LineWidth', 1);
ExactStress = plot(x, solu.sigma(x),'r-','LineWidth', 1);
legend([ApproxStress, ExactStress],'approx', 'exact')

[DNorm, ENorm] = calErrorNorm(Mesh, Line, E, d, solu);

disp(['DispNorm = ', num2str(DNorm, 20)]);
disp(['EnergyNorm = ', num2str(ENorm, 20)]);




