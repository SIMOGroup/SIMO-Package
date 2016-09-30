%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for one dimensional Euler-Bernoulli beam.
% (Rotation-free beam formulation due to high order continuity of 
% NURBS).
%

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
%
% Analytical Solution - Cantilever beam, point load
% --------------------------------------------------------------
% syms q E I L x F k
% uExact = F/(6*k)*(3*L*x.^2-x.^3);
% pretty(uExact);
%}
%{
% Analytical Solution - Cantilever beam, uniform load
% --------------------------------------------------------------
syms q E I L x F k
uExact = q/(24*k)*(x.^4 - 4*L*x.^3 + 6*L^2*x.^2);
pretty(uExact);
%}
% Material, geometry and force data

E = 100; % Young's modulus
h = 1; % thickness
%b = 1;
L = 8; % length of the beam
%I = b*h^3/12; % second moment of area
q = -1; % uniformly distributed loads
F = 10; % concentrated force
k = 1e6; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PRE-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control points

CtrlPts = zeros(4, 2);

CtrlPts(1 : 3, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2) = [L; 0; 0];

CtrlPts(4, :) = 1;

% knot vector

KntVect{1} = [0 0 1 1];
p=3;
nel=1;
kr=1;
Line = CreateNURBS(KntVect, CtrlPts);
Line = KRefine(Line, nel, p, p-kr);

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
df = @(x) q;
[KVals, FVals] = calcLocalStiffnessMatrices1DBeam(Mesh, Line, k, df);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals);
% assemble elementary stiffness matrix to the global matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

% Impose Dirichlet boundary conditions
disp([num2str(toc),'  Imposing Dirichlet boundary conditions'])
BCIdx = [1 2];
BCVal = [0 0]';     
FreeIdx = setdiff(1 : Mesh.NDof, BCIdx);
d = zeros(Mesh.NDof, 1);
d(BCIdx) = BCVal;
f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdx) * BCVal; 

% Impose Neumann boundary conditions
disp([num2str(toc),'  Imposing Neumann boundary conditions'])
f(end) = f(end) + F;

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
%plot(Line.CtrlPts3D(1, :), d, 'r.','MarkerSize', 20);
plot(Line.CtrlPts3D(1, :), d, 'ko-');
set(gcf, 'color', 'white');
% interpolate solution field
Pts{1} = linspace(0, 1, 10);
[C, U] = NURBSEval(Line, Pts, d);    
ApproxDispl = plot(C(1, :), U, 'bo-', 'LineWidth', 1);
%x  = C(1, :);
%solu = BeamExactSolution(L, k, q, F);
%ExactDispl = plot(x, solu.displ(x),'r*-', 'LineWidth', 1);
%legend([ApproxDispl, ExactDispl], 'IGA', 'Exact')
