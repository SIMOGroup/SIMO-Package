%---------------------------------------------------------------
% Insulated Wall Temperature
%---------------------------------------------------------------
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------PRE-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l1 = 0.75; %ft
l2 = 5/12; %ft
% control points

% line 1
CtrlPts{1} = zeros(4, 2);

CtrlPts{1}(1 : 3, 1) = [0; 0; 0];
CtrlPts{1}(1 : 3, 2) = [l1; 0; 0];

CtrlPts{1}(4, :) = 1;

% knot vector

KntVect1{1} = [0 0 1 1];

Line{1} = CreateNURBS(KntVect1, CtrlPts{1});

% line 2
CtrlPts{2} = zeros(4, 2);

CtrlPts{2}(1 : 3, 1) = [l1; 0; 0];
CtrlPts{2}(1 : 3, 2) = [l1 + l2; 0; 0];

CtrlPts{2}(4, :) = 1;

% knot vector

KntVect2{1} = [0 0 1 1];

Line{2} = CreateNURBS(KntVect2, CtrlPts{2});

Interfaces = struct;
Interfaces(1).Patch1 = 1;
Interfaces(1).Patch2 = 2;
Interfaces(1).Side1 = 2;
Interfaces(1).Side2 = 1;

Mesh = cell(1, numel(Line));
for iPtc = 1 : numel(Line)    
    Line{iPtc} = KRefine(Line{iPtc}, 5, 2, 0);
    Mesh{iPtc} = Mesh1D(Line{iPtc});
end

% figure
% hold on
% grid on
% axis equal
% daspect([1, 1, 1])

% for i = 1 : numel(Line)
%     PlotGeo(Line{i});
%     PlotKnts(Line{i});
%     PlotCtrlPts(Line{i});
%     PlotCtrlNet(Line{i});
% end

% Assemble siffnesss matrix
disp([num2str(toc),'  Assembling the system'])

[GNum, NDof] = MPInterfaceScalar1D(Interfaces, Mesh);

% -------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------PROCESSING--------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NCount = 0;

NPatch = numel(Mesh);

NEntsEmatK = zeros(1, NPatch);
NEntsGmatK = zeros(1, NPatch);

for iPtc = 1 : NPatch     
    NEntsEmatK(iPtc) = Mesh{iPtc}.NEN ^ 2;
    NEntsGmatK(iPtc) = NEntsEmatK(iPtc) * Mesh{iPtc}.NEl;
end
 
NEntsMult = sum(NEntsGmatK);

Rows = zeros(1, NEntsMult);
Cols = zeros(1, NEntsMult);
Vals = zeros(1, NEntsMult);
f = zeros(NDof, 1);

k(1) = 0.8; % Btu/hr-ft-°F
k(2) = 0.1; % Btu/hr-ft-°F
s{1} = @(x) 0;
s{2} = @(x) 0;

for iPtc = 1 : NPatch
    [KVals, FVals] = calcLocalConductionMatrices1D(Mesh{iPtc}, Line{iPtc},k(iPtc), s{iPtc}, 'PLANE WALL');
    [rsK, csK, vsK, f_locK] = convertToTripletStorage(Mesh{iPtc}, KVals, FVals);
    Rows(NCount + (1 : numel(rsK))) = GNum{iPtc}(rsK);
    Cols(NCount + (1 : numel(rsK))) = GNum{iPtc}(csK);
    Vals(NCount + (1 : numel(rsK))) = vsK;
    
    f(GNum{iPtc}) = f(GNum{iPtc}) + f_locK;
    
    NCount = NCount + numel(rsK);
end
K = sparse(Rows, Cols, Vals);

% Impose boundary conditions
disp([num2str(toc),'  Imposing boundary conditions'])

h(1) = 12; % Btu/hr-ft2-°F
h(2) = 2; % Btu/hr-ft2-°F

GDofs = GNum{2}(Mesh{2}.Boundary(2).Dofs);

K(1, 1) = K(1, 1) + h(1);
K(GDofs, GDofs) = K(GDofs, GDofs) + h(2);

T(1) = 3000; %°F
T(2) = 80; %°F

f(1) = h(1) * T(1);
f(GDofs) = h(2) * T(2);

% Solve the system
disp([num2str(toc),'  Solving the system'])
d = K \ f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------- POST-PROCESSING------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dptc = cell(1, numel(Mesh));
for i = 1 : numel(Mesh)
    dptc{i} = d(GNum{i});
end
% Plot control points
figure
hold on
grid on
set(gcf, 'color', 'white');
for i = 1 : numel(Mesh)
    plot(Line{i}.CtrlPts3D(1, :), dptc{i}, 'r.','MarkerSize',20);
    plot(Line{i}.CtrlPts3D(1, :), dptc{i}, 'k--');
end

% interpolate solution field
pts = linspace(0, 1, 1001);
C = cell(1, numel(Mesh));
U = cell(1, numel(Mesh));
for i = 1 : numel(Mesh)
    UXw = dptc{i}' .* Line{i}.CtrlPts4D(4, :);
    Uw = cat(1, Line{i}.CtrlPts4D, UXw);
    
    temp = BsplineEval(Line{i}.KntVect, Uw, {pts});
    Weights = temp(4, :);
    C{i} = temp(1, :) ./ Weights;
    U{i} = temp(5, :) ./ Weights;
end
 plot(cell2mat(C), cell2mat(U), 'k-', 'LineWidth', 1);

g1 = -k(1) * GradEval(Line{1}.KntVect, Line{1}.CtrlPts4D, Line{1}.uqKntVect, dptc{1}')
g2 = -k(2) * GradEval(Line{2}.KntVect, Line{2}.CtrlPts4D, Line{2}.uqKntVect, dptc{2}')


