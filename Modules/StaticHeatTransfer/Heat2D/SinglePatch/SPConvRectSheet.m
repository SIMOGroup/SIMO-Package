
% %
% % -------------------------------------------------------------                                 
% %    TWO-DIMENSIONAL HEAT TRANSFER WITH CONVECTION                              
% %                                                                               
% %    Test NAFEMS number T4                                                      
% %                                                                               
% %                                                                               
% %  description                                                                 
% %   -----------                                                                
% %                                                                              
% %             <- 0.6 m ->                                                      
% %            D     0°C   C  _____                                              
% %             .----------.       |                                             
% %            /|          |       |                                             
% %            /|          |       |                                             
% %            /|          |       |                                            
% %            /|          | 0°C   |                                            
% %            /|          |     1.0 m                                          
% %            /|          |       |                                             
% %            /|          |E__    |                                            
% %            /|          |  |    |                                           
% %            /|  100°C   | 0.2 m |                                            
% %             .----------. _|____|                                            
% %            A           B                                                    
% %                                                                             
% %     the internal heat production is zero                    
% %                                                              
% %     boundary conditions:                                                      
% %     ----------------------                                                                  
% %        - temperature imposed on the edge AB : T = 100°C                             
% %                                                                              
% %        - convection on edges BC and CD :                                     
% %            ambient temperature of 0°C                                        
% %                                                                              
% %        - edge AD is completely insulated : zero flux 
%         
% %     material properties:  
% %    -----------------------
% %       - isotropic material with thermal conductivity: 
% %           \kappa = 52 W/(°C m)
% %       - coefficient of heat transfer on edges BC, CD
% %           h = 750 W/(°C m^2)
% %                                                                               
% %   the known temperature at point E(0.6, 0.2) is 18.3°C                  
% %                                                              

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
b = 0.6; % m
h = 1; % m
h1 = 0.2; % m

CtrlPts = zeros(4, 2, 2);

CtrlPts(1 : 3, 1, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2, 1) = [b; 0; 0];

CtrlPts(1 : 3, 1, 2) = [0; h; 0];
CtrlPts(1 : 3, 2, 2) = [b; h; 0];

CtrlPts(4, :, :) = 1;

% knot vector

KntVect{1} = [0 0 1 1];
KntVect{2} = [0 0 1 1];

Surf = CreateNURBS(KntVect, CtrlPts);
Surf = KRefine(Surf, [1 1], [1 1], [0 0]);
% Surf = KRefine(Surf, [10 15], [4 4], [0 0]);

figure
hold on
% grid on
set(gcf,'color','white')
daspect([1 1 1])
axis([-0.1 b + 0.1 -0.1 h + 0.1]);
axis equal
axis off
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf);
PlotCtrlNet(Surf);

% --------------------------------------------------------------
Mesh = Mesh2D(Surf, 'ScalarField');
% material properties
ka = 52; % W/(m * ^\circ C)
s = @(x, y) 0;

[KVals, FVals4] = calcLocalConductionMatrices2D(Mesh, Surf, ka, s);
[Rows, Cols, Vals, f] = convertToTripletStorage(Mesh, KVals, FVals4);

% Impose convection boundary condition
disp([num2str(toc),'  Imposing convection boundary condition'])
hCoeff = 750; 
Ta = 0; 
Ref1 = 4; % boundary reference
[HVals4, FVals4] = calcLocalHeatTransferMatrices2D(Surf, Mesh, hCoeff, Ta, Ref1);
[RowsH4, ColsH4, ValsH4, f1] = convertToTripletStorage(Mesh.Boundary(Ref1), HVals4, FVals4);

Ref2 = 2;
[HVals2, FVals2] = calcLocalHeatTransferMatrices2D(Surf, Mesh, hCoeff, Ta, Ref2);
[RowsH2, ColsH2, ValsH2, f2] = convertToTripletStorage(Mesh.Boundary(Ref2), HVals2, FVals2);

Rs = [Rows; Mesh.Boundary(Ref1).Dofs(RowsH4); Mesh.Boundary(Ref2).Dofs(RowsH2)];
Cs = [Cols; Mesh.Boundary(Ref1).Dofs(ColsH4); Mesh.Boundary(Ref2).Dofs(ColsH2)];
Vs = [Vals; ValsH4; ValsH2];

% Convert triplet data to sparse matrix
K = sparse(Rs, Cs, Vs);
clear Rs Cs Vs

f(Mesh.Boundary(4).Dofs) = f(Mesh.Boundary(4).Dofs) + f1;
f(Mesh.Boundary(2).Dofs) = f(Mesh.Boundary(2).Dofs) + f2;

% Impose essential boundary conditions
disp([num2str(toc),'  Imposing essential boundary conditions'])
TAB = 100; %°C
BCIdx = Mesh.Boundary(3).Dofs;
BCVal = repmat(TAB, numel(BCIdx), 1);

FreeIdx = setdiff(1 : Mesh.NDof, BCIdx);

d = zeros(Mesh.NDof, 1);
d(BCIdx) = BCVal;

f(FreeIdx) = f(FreeIdx) - K(FreeIdx, BCIdx) * BCVal; 

% Solve the system
disp([num2str(toc),'  Solving the system'])
d(FreeIdx) = K(FreeIdx, FreeIdx) \ f(FreeIdx);

% Plot heat flux at gauss points
disp([num2str(toc),'  Ploting heat flux'])
PlotFlux2D(Surf, Mesh, d, ka)

% Export result to *.vtk format
disp([num2str(toc),'  Exporting result to *.vtk format'])
ParaPts = {linspace(0, 1, 101), linspace(0, 1, 101)};

SPToVTK(Surf, d, ParaPts, 'SPConvRectSheet', 'Temperature')

[C, F] = NURBSEval(Surf, {1, 0.2}, d);
sprintf('benchmark = %g', F)

