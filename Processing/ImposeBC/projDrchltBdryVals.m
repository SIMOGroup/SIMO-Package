function [Coeffs, GDofs] = projDrchltBdryVals(NURBS, Mesh, h, Refs, LAB, varargin)
% function [Coeffs, GDofs] = projDrchltBdryVals(NURBS, Mesh, h, Refs, LAB, varargin)
% Evaluate coefficent values for imposing Dirichlet boundary condition of
% 2D and 3D problem (project Dirichlet Boundary Values)
% -------------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure (single patch) or cell of NURBS structures
%       (multiple patches)
%       Mesh: Mesh structure (single patch) or cell of Mesh structures
%       (multiple patches)
%       h: boundary function, Ex: h = @(x, y) a * x + b * y
%           h = @(x, y) 0 correspond to homogeneous Dirichlet B.C
%       Refs: referenced indices to indicate boundary curves,
%       Refs = [Ref_1, Ref_2,...,Ref_n],
%       ******************************************************
%       *                                                    *
%       *                   Ref = 4                          *
%       *                   (v = 1)                          *
%       *               o------------o              v        *
%       *               |            |              ^        *
%       *               |            |              |        *
%       *      Ref = 1  |            | Ref = 2      |        *
%       *       (u = 0) |            |   (u = 1)    +----> u *
%       *               |            |                       *
%       *               o------------o                       *
%       *                   Ref = 3                          *
%       *                   (v = 0)                          *
%       ******************************************************
%       LAB: 'TEMP' for thermal problem
%            'UX', 'UY' or 'UZ' for structural problem
%            'PLATE' for plate problem
%       varargin (optional):
%       - if varargin = GNum: single patch problem
%       with coupled dofs.
%       - if varargin = {GNum, Boundaries}: multiple patches problem.
% -------------------------------------------------------------------------
% Output:
%       GDofs: indices of degrees of freedom of the given boundary
%       Coeffs: evaluated coefficent values corresponding to these degrees
%       of freedom
% -------------------------------------------------------------------------

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

GDofs = []; GARows = []; GACols = []; GAVals = []; GFVals = []; GFIdcs = [];
assert(isa(h, 'function_handle'), 'h must be a function handle')
assert(ischar(LAB), 'You must specify the LABEL keyword')
for iRef = 1 : numel(Refs) % Loop over boundaries
    if nargin == 7 % multiple patches problem (NURBS, Mesh, h, Refs, LAB, GNum, Boundaries)
        GNum = varargin{1};
        Boundaries = varargin{2};
        NPatches = numel(Boundaries(Refs(iRef)).Patches);
        for BdrySide = 1 : NPatches % Loop over patches of the boundary
            iPtc = Boundaries(Refs(iRef)).Patches(BdrySide);
            iSide = Boundaries(Refs(iRef)).Sides(BdrySide);
            [GARows, GACols, GAVals, GFIdcs, GFVals, GDofs] = getDrchltBdryData(NURBS{iPtc}, Mesh{iPtc}, iSide, h, LAB, GARows, GACols, GAVals, GFIdcs, GFVals, GDofs, GNum{iPtc});
        end
    else % single patch problem
        iSide = Refs(iRef);
        [GARows, GACols, GAVals, GFIdcs, GFVals, GDofs] = getDrchltBdryData(NURBS, Mesh, iSide, h, LAB, GARows, GACols, GAVals, GFIdcs, GFVals, GDofs, varargin{:});
    end
end
% Convert triplet data to sparse matrix
% --------------------------------------------------------------------
A = sparse(GARows, GACols, GAVals);
F = sparse(GFIdcs, ones(numel(GFIdcs), 1), GFVals);
% --------------------------------------------------------------------
Coeffs = full(A(GDofs, GDofs) \ F(GDofs));
end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Auxiliary functions
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [GARows, GACols, GAVals, GFIdcs, GFVals, GDofs] =...
    getDrchltBdryData(NURBS, Mesh, iSide, h, LAB, GARows, GACols, GAVals, GFIdcs, GFVals, GDofs, varargin)
MeshBdry = Mesh.Boundary(iSide);
NURBSBdry = NURBSBoundary(NURBS, iSide);
if strcmp(LAB, 'TEMP') || strcmp(LAB, 'PLATE') % thermal or plate problem
    if nargin == 12 && isvector(varargin{:}) % coupling dofs or multiple patches of thermal or plate problem
        LDofs = varargin{:}(MeshBdry.Dofs);
    else
        LDofs = MeshBdry.Dofs;
    end
elseif strcmp(LAB, 'UX')% structural problem
    if nargin == 12 && isvector(varargin{:}) % coupling dofs or multiple patches of structural problem
        LDofs = varargin{:}(MeshBdry.CompDofs{1});
    else
        LDofs = MeshBdry.CompDofs{1};
    end
elseif strcmp(LAB, 'UY')
    if nargin == 12 && isvector(varargin{:})
        LDofs = varargin{:}(MeshBdry.CompDofs{2});
    else
        LDofs = MeshBdry.CompDofs{2};
    end
elseif strcmp(LAB, 'UZ')
    if nargin == 12 && isvector(varargin{:})
        LDofs = varargin{:}(MeshBdry.CompDofs{3});
    else
        LDofs = MeshBdry.CompDofs{3};
    end
end
if NURBSBdry.Dim == 1
    [LAVals, LFVals] = applyL2Proj2D(NURBSBdry, MeshBdry, h);
elseif NURBSBdry.Dim ==2
    [LAVals, LFVals] = applyL2Proj3D(NURBSBdry, MeshBdry, h);
end
J = repmat(1 : MeshBdry.NEN, MeshBdry.NEN, 1);
I = J';
ii = MeshBdry.El(:, I(:))';
jj = MeshBdry.El(:, J(:))';
LARows = LDofs(ii(:));
LACols = LDofs(jj(:));
% ---------------------------------------------
tmp = MeshBdry.El'; % indices of local force vector
LFIdcs = LDofs(tmp(:));

GARows = cat(1, GARows, LARows);
GACols = cat(1, GACols, LACols);
GAVals = cat(1, GAVals, LAVals);
GFIdcs = cat(1, GFIdcs, LFIdcs);
GFVals = cat(1, GFVals, LFVals);
GDofs = union(GDofs, LDofs);
end
% -----------------------------------------------------------------------
% L2 Projection: project an arbitrary function h \in L2(\Gamma) into a
% finite element space Vh \in L2(\Gamma), \Gamma is the boundary of the domain
% Input:
%       NURBS: NURBS structure of the boundary curve
%       Mesh: Mesh structure of the boundary curve
%       h: boundary function, Ex: h = @(x, y) a * x + b * y
%           h = @(x, y) 0 correspond to homogeneous Dirichlet B.C
% Output:
%       LAVals: local mass matrices
%       LFVals: local force vectors
function [LAVals, LFVals] = applyL2Proj2D(NURBS, Mesh, h)
LAVals = [];
LFVals = [];
[CtrlPts, Weights] = convertTo2DArrays(NURBS);
NGPs = NURBS.Order + 1;
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order, NURBS.NCtrlPts, NURBS.KntVect{1}, 1, NGPs, Mesh.NEl);
for e = 1 : Mesh.NElDir(1)
    LA = zeros(Mesh.NEN);
    LF = zeros(Mesh.NEN, 1);
    for qx = 1 : NGPs
        N0 = Nx(e, qx, :, 1);
        N1 = Nx(e, qx, :, 2);
        [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0(:)', N1(:)');
        % gradient of mapping from parameter space to physical space
        dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
        % compute the jacobian of physical and parameter domain mapping
        J1 = norm(dxdxi);
        LA = LA + R0' * R0 * J1 * Jx(e) * Wx(qx);
        
        Pts = R0 * CtrlPts(Mesh.El(e, :), :);
        x = Pts(1); y = Pts(2);
        LF = LF + h(x, y) * R0' * J1 * Jx(e) * Wx(qx);
    end
    LAVals = cat(1, LAVals, LA(:));
    LFVals = cat(1, LFVals, LF);
end
end
% ------------------------------------------------------------------------
function [LAVals, LFVals] = applyL2Proj3D(NURBS, Mesh, h)
LAVals = [];
LFVals = [];
[CtrlPts, Weights] = convertTo2DArrays(NURBS);
NGPs = NURBS.Order + 1;
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));
N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);
for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        e = sub2ind(Mesh.NElDir, ex, ey);
        LA = zeros(Mesh.NEN);
        LF = zeros(Mesh.NEN, 1);
        for qy = 1 : NGPs(2)
            for qx = 1 : NGPs(1)
                k = 1;
                for j = 1 : NURBS.Order(2) + 1
                    for i = 1 : NURBS.Order(1) + 1
                        N0(k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1);
                        N1(1, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 1);
                        N1(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 2);
                        k = k + 1;
                    end
                end
                [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0, N1);
                
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % Gradient of mapping from parameter space to physical
                % space
                dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
                % compute the jacobian of physical and parameter domain mapping
                t1 = dxdxi(1, :);
                t2 = dxdxi(2, :);
                n = cross(t1, t2);
                J1 = norm(n);
                LA = LA + R0' * R0 * J1 * J2 * W;
                Pts = R0 * CtrlPts(Mesh.El(e, :), :);
                x = Pts(1); y = Pts(2);
                LF = LF + h(x, y) * R0' * J1 * J2 * W;
            end
        end
        LAVals = cat(1, LAVals, LA(:));
        LFVals = cat(1, LFVals, LF);
    end
end
end