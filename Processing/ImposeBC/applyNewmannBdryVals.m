function [Vals, GDofs] = applyNewmannBdryVals(NURBS, Mesh, g, Refs, LAB, varargin)
% [Vals, GDofs] = applyNewmannBdryVals(NURBS, Mesh, g, Refs, LAB, varargin)
% Evaluate values for applying Newmann boundary condition
% ------------------------------------------------------------
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
%       *                  (eta = 1)                         *
%       *            P3 o------------o P4           v        *
%       *               |            |              ^        *
%       *               |            |              |        *
%       *      Ref = 1  |            | Ref = 2      |        *
%       *      (xi = 0) |            | (xi = 1)     +----> u *
%       *               |            |                       *
%       *            P1 o------------o P2                    *
%       *                   Ref = 3                          *
%       *                  (eta = 0)                         *
%       ******************************************************
%       LAB: 'HFLUX' Heat Fluxes for thermal problem
%            'FX', 'FY' or 'FZ' for structural problem
%            'PRESS' Pressure for structural problem
%       varargin (optional):
%       - if varargin = GNum: single patch problem
%       with coupled dofs.
%       - if varargin = {GNum, Boundaries}: multiple patches problem.
% ------------------------------------------------------------------
% Output:
%       GDofs: indices of degrees of freedom of the given boundary
%       Vals: evaluated coefficent values corresponding to these degrees
%       of freedom
% ------------------------------------------------------------------

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

GDofs = []; GVals = []; GIdcs = [];
assert(isa(g, 'function_handle'), 'g must be a function handle')
assert(ischar(LAB), 'You must specify the LABEL keyword')
for iRef = 1 : numel(Refs) % Loop over boundaries
    if nargin == 7 % multiple patches problem (NURBS, Mesh, h, Refs, LAB, GNum, Boundaries)
        GNum = varargin{1};
        Boundaries = varargin{2};
        NPatches = numel(Boundaries(Refs(iRef)).Patches);
        for BdrySide = 1 : NPatches % Loop over patches of the boundary
            iPtc = Boundaries(Refs(iRef)).Patches(BdrySide);
            iSide = Boundaries(Refs(iRef)).Sides(BdrySide);
            [GIdcs, GVals, GDofs] = getNewmannBdryData(NURBS{iPtc}, Mesh{iPtc}, iSide, g, LAB, GIdcs, GVals, GDofs, GNum{iPtc});
        end
    else % single patch problem
        iSide = Refs(iRef);
        [GIdcs, GVals, GDofs] = getNewmannBdryData(NURBS, Mesh, iSide, g, LAB, GIdcs, GVals, GDofs, varargin{:});
    end
end
Vals = sparse(GIdcs, ones(numel(GIdcs), 1), GVals);
Vals = full(Vals(GDofs));
end
function [GIdcs, GVals, GDofs] = getNewmannBdryData(NURBS, Mesh, iSide, g, LAB, GIdcs, GVals, GDofs, varargin)
% function [GVals, GDofs] = getNewmannBdryData(NURBS, Mesh, iSide, g, LAB, GVals, GDofs, varargin)
MeshBdry = Mesh.Boundary(iSide);
NURBSBdry = NURBSBoundary(NURBS, iSide);
if strcmp(LAB, 'HFLUX') || strcmp(LAB, 'hflux') || strcmp(LAB, 'HFlux') || strcmp(LAB, 'PRES') || strcmp(LAB, 'pres')
    if nargin == 9 && isvector(varargin{:}) % coupling dofs or multiple patches of thermal or plate problem
        LDofs = varargin{:}(MeshBdry.Dofs);
    else
        LDofs = MeshBdry.Dofs;
    end
elseif strcmp(LAB, 'FX') || strcmp(LAB, 'fx') || strcmp(LAB, 'Fx')
    if nargin == 9 && isvector(varargin{:}) % coupling dofs or multiple patches of structural problem
        LDofs = varargin{:}(MeshBdry.CompDofs{1});
    else
        LDofs = MeshBdry.CompDofs{1};
    end
elseif strcmp(LAB, 'FY') || strcmp(LAB, 'fy') || strcmp(LAB, 'Fy')
    if nargin == 9 && isvector(varargin{:})
        LDofs = varargin{:}(MeshBdry.CompDofs{2});
    else
        LDofs = MeshBdry.CompDofs{2};
    end
elseif strcmp(LAB, 'FZ') || strcmp(LAB, 'fz') || strcmp(LAB, 'Fz')
    if nargin == 9 && isvector(varargin{:})
        LDofs = varargin{:}(MeshBdry.CompDofs{3});
    else
        LDofs = MeshBdry.CompDofs{3};
    end
end
if NURBSBdry.Dim == 1
    if strcmp(LAB, 'PRES')
        LVals = calcPressureVals2D(NURBSBdry, MeshBdry, g, iSide);
        tmp = [MeshBdry.El, MeshBdry.El + NURBSBdry.NNP]';
    else
        LVals = calcNewmannVals2D(NURBSBdry, MeshBdry, g);
        tmp = MeshBdry.El';
    end
elseif NURBSBdry.Dim ==2
    if strcmp(LAB, 'PRES')
        LVals = calcPressureVals3D(NURBSBdry, MeshBdry, g);
        tmp = [MeshBdry.El, MeshBdry.El + NURBSBdry.NNP, MeshBdry.El + NURBSBdry.NNP * 2]';
    else
        LVals = calcNewmannVals3D(NURBSBdry, MeshBdry, g);
        tmp = MeshBdry.El';
    end
end
assert(iscolumn(LDofs), 'Data must be stored in column-vector format')
LIdcs = LDofs(tmp(:));
GIdcs = cat(1, GIdcs, LIdcs);
GVals = cat(1, GVals, LVals);
GDofs = union(GDofs, LDofs);
end
function LVals = calcNewmannVals2D(NURBS, Mesh, g)
% function LVals = calcNewmannVals2D(NURBS, Mesh, g)
NGPs = NURBS.Order + 1;
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order, NURBS.NCtrlPts, NURBS.KntVect{1}, 1, NGPs, Mesh.NEl);
LFVals = zeros(Mesh.NEN, Mesh.NEl);
[CtrlPts, Weights] = convertTo2DArrays(NURBS);
for e = 1 : Mesh.NEl
    LF = zeros(Mesh.NEN, 1);
    for qx = 1 : NGPs
        N0 = Nx(e, qx, :, 1);
        N1 = Nx(e, qx, :, 2);
        [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0(:)', N1(:)');
        % gradient of mapping from parameter space to
        % physical space
        dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
        % compute the jacobian of physical and parameter
        % domain mapping
        J1 = norm(dxdxi);
        Pts = R0 * CtrlPts(Mesh.El(e, :), :);
        x = Pts(1); y = Pts(2);
        LF = LF + Wx(qx) * R0' * g(x, y) * J1 * Jx(e);
    end
    LFVals(:, e) = LF;
end
LVals = LFVals(:);
end
function LVals = calcNewmannVals3D(NURBS, Mesh, g)
% function Vals = calcNewmannVals3D(NURBS, Mesh, g)
NGPs = NURBS.Order + 1;
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));
[CtrlPts, Weights] = convertTo2DArrays(NURBS);
LFVals = zeros(Mesh.NEN, Mesh.NEl);
for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        e = sub2ind(Mesh.NElDir, ex, ey);
        LF = zeros(Mesh.NEN, 1);
        for qy = 1 : NGPs(2)
            for qx = 1 : NGPs(1)
                N0x = Nx(ex, qx, :, 1);
                N1x = Nx(ex, qx, :, 2);
                N0y = Ny(ey, qy, :, 1);
                N1y = Ny(ey, qy, :, 2);
                N0 = bsxfun(@times, N0x(:), N0y(:)');
                N11 = bsxfun(@times, N1x(:), N0y(:)');
                N12 = bsxfun(@times, N0x(:), N1y(:)');
                [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0(:)', [N11(:)'; N12(:)']);
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % gradient of mapping from parameter space to physical space
                dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
                % compute the jacobian of physical and parameter domain mapping
                t1 = dxdxi(1, :);
                t2 = dxdxi(2, :);
                n = cross(t1, t2);
                J1 = norm(n);
                Pts = R0 * CtrlPts(Mesh.El(e, :), :);
                x = Pts(1); y = Pts(2); z = Pts(3);
                LF = LF + R0' * g(x, y, z) * J1 * J2 * W;
            end
        end
        LFVals(:, e) = LF;
    end
end
LVals = LFVals(:);
end
function LVals = calcPressureVals2D(NURBS, Mesh, g, iSide)
% function Vals = calcPressureVals2D(NURBS, Mesh, g, iSide)
NGPs = NURBS.Order + 1; % number of gauss points = p + 1
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order, NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs, Mesh.NElDir(1));
[CtrlPts, Weights] = convertTo2DArrays(NURBS);
LFVals = zeros(Mesh.NEN * 2, Mesh.NEl);
for e = 1 : Mesh.NElDir(1)
    LF = zeros(Mesh.NEN * 2, 1);
    for qx = 1 : NGPs
        N0 = Nx(e, qx, :, 1);
        N1 = Nx(e, qx, :, 2);
        [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0(:)', N1(:)');
        % gradient of mapping from parameter space to physical space
        dxdxi = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
        % compute the jacobian of physical and parameter domain mapping
        J1 = norm(dxdxi);
        axis = floor((iSide + 1) / 2);
        n = [0, -1; 1, 0] * dxdxi'; % dxdxi: tangential vector
        n = n / norm(n);
        if axis == 1
            n = -n;
        end
        if rem(iSide, 2) ~= 0 % side 1, 3
            n = -n;
        end
        Pts = R0 * CtrlPts(Mesh.El(e, :), :);
        x = Pts(1); y = Pts(2);
        % t = dxdxi / norm(dxdxi);
        % quiver(x, y, t(1), t(2), 'k');
        % quiver(x, y, n(1), n(2), 'k');
        p = repmat(g(x, y), 2, 1) .* n;
        quiver(x, y, p(1), p(2), 0.02, 'k');
        R0Mat = zeros(2, Mesh.NEN * 2); % matrix of basis functions
        R0Mat(1, 1 : Mesh.NEN) = R0;
        R0Mat(2, Mesh.NEN + 1 : 2 * Mesh.NEN) = R0;
        LF = LF + Wx(qx) * R0Mat' * p * J1 * Jx(e);
    end
    LFVals(:, e) = LF;
end
LVals = LFVals(:);
end
function LVals = calcPressureVals3D(NURBS, Mesh, g)
% function Vals = calcPressureVals3D(NURBS, Mesh, g)
NGPs = NURBS.Order + 1; % number of gauss points = p + 1
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 1, NGPs(2), Mesh.NElDir(2));
[CtrlPts, Weights] = convertTo2DArrays(NURBS);
LFVals = zeros(Mesh.NEN * 3, Mesh.NEl);
for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        LF = zeros(Mesh.NEN * 3, 1);
        e = sub2ind(Mesh.NElDir, ex, ey);
        for qy = 1 : NGPs(2)
            for qx = 1 : NGPs(1)
                N0x = Nx(ex, qx, :, 1);
                N1x = Nx(ex, qx, :, 2);
                N0y = Ny(ey, qy, :, 1);
                N1y = Ny(ey, qy, :, 2);
                N0 = bsxfun(@times, N0x(:), N0y(:)');
                N11 = bsxfun(@times, N1x(:), N0y(:)');
                N12 = bsxfun(@times, N0x(:), N1y(:)');
                [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0(:)', [N11(:)'; N12(:)']);
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % gradient of mapping from parameter space to physical space
                dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
                % compute the jacobian of physical and parameter domain mapping
                t1 = dxdxi(1, :);
                t2 = dxdxi(2, :);
                n = cross(t1, t2);
                J1 = norm(n);
                n = n / J1;
                Pts = R0 * CtrlPts(Mesh.El(e, :), :);
                x = Pts(1); y = Pts(2); z = Pts(3);
                %quiver3(x, y, z, t1(1), t1(2), t1(3))
                %quiver3(x, y, z, t2(1), t2(2), t2(3))
                %quiver3(x, y, z, n(1), n(2), n(3))
                p = repmat(g(x, y, z), 3, 1) .* n';
                quiver3(x, y, z, p(1), p(2), p(3))
                R0Mat = zeros(3, 3 * Mesh.NEN); % matrix of basis functions
                R0Mat(1, 1 : Mesh.NEN) = R0;
                R0Mat(2, Mesh.NEN + 1 : 2 * Mesh.NEN) = R0;
                R0Mat(3, 2 * Mesh.NEN + 1 : 3 * Mesh.NEN) = R0;
                LF = LF + W * R0Mat' * p * J1 * J2;
            end
        end
        LFVals(:, e) = LF;
    end
end
LVals = LFVals(:);
end