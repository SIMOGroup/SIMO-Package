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
E = 1; % Young's Modulus

% physical parameters
A = 1; % Area of the bar
L = 1; % Length of the bar

% essential boundary condition
u0 = 0;
u1 = 1;

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

Line = HRefine(Line, 1, [0.42, 0.5, 0.58]);
% % line = HRefine(line, 2, 0.5);
% line = KRefine(line, 2, 4, 2);
Line = KRefine(Line, 3, 3, 1);

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
disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])
alfa = 50;
bfx = @(x) (x >= 0.42 && x <= 0.58) *...
    (2 * alfa ^ 2 - 4 * (alfa ^ 2 * (x - 0.5)) ^ 2) *...
    exp(-(alfa * (x - 0.5)) ^ 2);

[KVals, FVals] = calcLocalStiffnessMatrices1D(Mesh, Line, E, A, bfx);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals);

% assemble elementary stiffness matrix to the global matrix
K = sparse(Rows, Cols, Vals);

% Impose Dirichlet boundary conditions
disp([num2str(toc),'  Imposing Dirichlet boundary conditions'])
BCIdx = [1 Line.NNP];
BCVal = [u0; u1];
FreeIdx = setdiff(1 : Mesh.NDof, BCIdx);
d = zeros(Mesh.NDof, 1);
d(BCIdx) = BCVal;
f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdx) * BCVal; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdx) = K(FreeIdx, FreeIdx) \ f(FreeIdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------- POST-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% interpolate solution field
Pts = linspace(0, 1, 1001);

UXw = d' .* Line.CtrlPts4D(4, :);
Uw = cat(1, Line.CtrlPts4D, UXw);

temp = BsplineEval(Line.KntVect, Uw, {Pts});
Weights = temp(4, :);
C = temp(1, :) ./ Weights;
U = temp(5, :) ./ Weights;

figure    
plot(C, U, 'k-', 'LineWidth', 1);

% plot exact solution of displacement
% % x  = pts;
uExact = @(x) x + exp(-(alfa * (x - 0.5)) .^ 2);
plot(Pts, uExact(Pts),'r-', 'LineWidth',1)

figure
hold on
set(gcf, 'color', 'white');
% evaluate strain
g = GradEval(Line.KntVect, Line.CtrlPts4D, {Pts}, d');

plot(C, g * E, 'k-','LineWidth',1);

% plot exact solution of stress

sigmaExact = @(x) 1 - 2 .* alfa ^ 2 .* (x - 0.5) .* exp(-(alfa .* (x - 0.5)) .^ 2);
plot(Pts, sigmaExact(Pts),'r-','LineWidth',1)

