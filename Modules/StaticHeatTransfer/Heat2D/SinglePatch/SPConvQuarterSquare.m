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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control points
tic;
r = 0.4;
a = 1;

CtrlPts1 = zeros(4, 4);
CtrlPts1(1 : 3, 1) = [-a; a; 0];
CtrlPts1(1 : 3, 2) = [-a; 0; 0];
CtrlPts1(1 : 3, 3) = [-a; 0; 0];
CtrlPts1(1 : 3, 4) = [0; 0; 0];

CtrlPts1(4, :) = 1;

Curv2 = NURBSCirc((r + a) / 2, [0, a, 0], pi, 3/2*pi);
Curv2 = HRefine(Curv2, 1, 1/2);
CtrlPts2 = Curv2.CtrlPts4D;
clear Curv2

Curv3 = NURBSCirc(r, [0, a, 0], pi, 3/2*pi);
Curv3 = HRefine(Curv3, 1, 1/2);
CtrlPts3 = Curv3.CtrlPts4D;
clear Curv3

KntVect{1} = [0 0 0 1/2 1 1 1];
KntVect{2} = [0 0 0 1 1 1];

CtrlPtsSurf = cat(3, CtrlPts1, CtrlPts2, CtrlPts3);

Surf = CreateNURBS(KntVect, CtrlPtsSurf);
Surf = KRefine(Surf, [5 5], [2 2], [1 1]);

figure
hold on
grid on
set(gcf,'color','white')
daspect([1 1 1])
axis equal
% axis off

PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf);
PlotCtrlNet(Surf);

Mesh = Mesh2D(Surf, 'ScalarField');
% material properties
ka = 1; % W/(m * ^\circ C)
s = @(x, y) -4;

% Assemble conduction matrix
disp([num2str(toc),'  Assembling the system'])
[KVals, FVals] = calcLocalConductionMatrices2D(Mesh, Surf, ka, s);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals);

% Impose convection boundary condition
disp([num2str(toc),'  Imposing convection boundary condition'])
hCoeff = 3; % W/m^2 \circ C
Ta = -0.32/hCoeff; % \circ C

Ref1 = 4;
[HVals, FVals] = calcLocalHeatTransferMatrices2D(Surf, Mesh, hCoeff, Ta, Ref1);
[RowsH, ColsH, ValsH, f1] = convertToTripletStorage(Mesh.Boundary(Ref1), HVals, FVals);

Rs = [Rows; Mesh.Boundary(Ref1).Dofs(RowsH)];
Cs = [Cols; Mesh.Boundary(Ref1).Dofs(ColsH)];
Vs = [Vals; ValsH];

f(Mesh.Boundary(4).Dofs) = f(Mesh.Boundary(4).Dofs) + f1;

% Convert triplet data to sparse matrix
K = sparse(Rs, Cs, Vs);
clear Rs Cs Vs f1

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary condition'])
T0 = @(x, y) (x ^ 2 + (y - 1) ^ 2);

[BCVals, BCIdcs] = projDrchltBdryVals(Surf, Mesh, T0, 3, 'TEMP');

FreeIdx = setdiff(1 : Mesh.NDof, BCIdcs);

d = zeros(Mesh.NDof, 1);
d(BCIdcs) = BCVals;

f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdcs) * BCVals; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdx) = K(FreeIdx, FreeIdx) \ f(FreeIdx);

% % Plot heat flux at gauss points
% disp([num2str(toc),'  Ploting heat flux'])
% PlotFlux2D(Surf, Mesh, d, ka)

% Evaluate temperture along diagonal line
[C, F] = NURBSEval(Surf, {0.5, linspace(0, 1, 1001)}, d);
alfa = atan2(C(2, :) - 1, C(1, :));
plot(C(1, :), C(2, :), 'LineWidth', 1.5);
figure
hold on
x = C(1, :);
y = C(2, :);
plot(x, x.^2 + (y - 1).^2, 'r', 'LineWidth', 1.5);
plot(x, F, 'b');
legend('exact', 'approx')

% Export result to *.vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101)};
filename = 'SPConvQuarterSquare';
fieldname = 'Temperature';
SPToVTK(Surf, d, ParaPts, filename, fieldname)

