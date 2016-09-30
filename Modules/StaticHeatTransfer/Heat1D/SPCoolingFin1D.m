% cooling fin
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

clear
close all
clc
tic
%-----------------------------------------
%  input data for coefficients of the ODE
%-----------------------------------------
ka = 385;
a = 100;
t = 1;
L = 5;
h = 25;
T0 = 100;
Ta = 20;
P = 2 * (L + t);
A = L * t;

Line = NURBSLine([0, 0, 0], [a, 0, 0]);

Line = HRefine(Line, 1, 0.2);
Line = KRefine(Line, 3, 4, 0);

figure
hold on
axis equal
daspect([1, 1, 1])
PlotGeo(Line);
PlotKnts(Line);
PlotCtrlPts(Line);
PlotCtrlNet(Line);

Mesh = Mesh1D(Line);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------PROCESSING--------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'  Assembling the system'])
s = @(x) 0; % heat source

[KVals, FVals] = calcLocalConductionMatrices1D(Mesh, Line, ka, s, 'FIN', P, A, h, Ta);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals);

% assemble elementary stiffness matrix to the global matrix
K = sparse(Rows, Cols, Vals);
clear Rows Cols Vals

% introducing surface heat transfer matrix at the right end
K(end, end) = K(end, end) + A * h;

% appling ambient-temperature load
f(end) = f(end) + A * h * Ta;
% Impose Dirichlet boundary conditions
disp([num2str(toc),'  Imposing Dirichlet boundary conditions'])
BCIdx = 1;
BCVal = T0;
FreeIdx = setdiff(1 : Mesh.NDof, BCIdx);
d = zeros(Mesh.NDof, 1);
d(BCIdx) = BCVal;
f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdx) * BCVal; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdx) = K(FreeIdx, FreeIdx) \ f(FreeIdx);

figure
hold on

% Evaluate the exact solution
pts2 = linspace(0, a, 1000);
xi = pts2 / a;
m = sqrt((h * P * a ^ 2) / (ka * A));
N = h * a / ka;
theta = (m * cosh(m * (1 - xi)) + N * sinh(m * (1 - xi))) /...
    (m * cosh(m) + N * sinh(m));
TEx = theta*(T0 - Ta) + Ta;

Pts{1} = linspace(0, 1, 1001);
[C, ApproxTemp] = NURBSEval(Line, Pts, d);    

plot1(1) = plot(C(1, :), ApproxTemp, 'k');
plot1(2) = plot(pts2, TEx,'k');
set(plot1(1),'DisplayName','Approx','LineStyle','-.');
set(plot1(2),'DisplayName','Exact');
legend('Approx','Exact')













