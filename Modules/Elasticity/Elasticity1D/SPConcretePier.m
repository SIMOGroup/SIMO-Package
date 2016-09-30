% Concrete Pier
% reference: Chapter 7, Reddy, Principles of Continuum Mechanics, Cambridge
% University Press.

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
E = @(x) 28e6; % Young's Modulus, KN/m^2

% physical parameters
A = @(x) 1 /4 * (1 + x); % Area of the bar, m^2
L = 2; % Length of the bar

% essential boundary condition
u0 = 0;

% natural boundary conditions
P = 20 * 0.25; % KN
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

Line = CreateNURBS(KntVect, CtrlPts);
Line = KRefine(Line, 2, 2, 1);

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
df = @(x) 6.25 * (1 + x); % body force per unit length

[KVals, FVals] = calcLocalStiffnessMatrices1D(Mesh, Line, E, A, df);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals);

% assemble elementary stiffness matrix to the global matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

% Impose Dirichlet boundary conditions
disp([num2str(toc),'  Imposing Dirichlet boundary conditions'])
BCIdx = Mesh.NDof;
BCVal = u0;
FreeIdx = setdiff(1 : Mesh.NDof, BCIdx);
d = zeros(Mesh.NDof, 1);
d(BCIdx) = BCVal;
f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdx) * BCVal; 

% Impose Neumann boundary conditions
disp([num2str(toc),'  Imposing Neumann boundary conditions'])
f(1) = f(1) + P;

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
% plot(Line.CtrlPts3D(1, :), d, 'r.','MarkerSize', 20);
% plot(Line.CtrlPts3D(1, :), d, 'k--');
set(gcf, 'color', 'white');
% interpolate solution field
Pts{1} = linspace(0, 1, 1001);
[C, U] = NURBSEval(Line, Pts, d);    
ApproxDispl = plot(C(1, :), U, 'k-', 'LineWidth', 1);
x = C(1, :);
% Evaluate exact solution
solu = SPConcretePierExactSolution(E);
ExactDispl = plot(x, solu.displ(x),'r-', 'LineWidth', 1);
legend([ApproxDispl, ExactDispl], 'approx', 'exact')

figure
hold on
title('Stress')
set(gcf, 'color', 'white');
g = GradEval(Line.KntVect, Line.CtrlPts4D, Pts, d'); % evaluate strain
ApproxStress = plot(x, g * E(x), 'k-','LineWidth', 1);
ExactStress = plot(x, solu.sigma(x),'r-','LineWidth', 1);
legend([ApproxStress, ExactStress],'approx', 'exact')

[DNorm, ENorm] = calErrorNorm(Mesh, Line, E, d, solu);
disp(['DispNorm = ', num2str(DNorm, 20)]);
disp(['EnergyNorm = ', num2str(ENorm, 20)]);
