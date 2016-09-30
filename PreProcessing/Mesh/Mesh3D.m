function Mesh = Mesh3D(NURBS, Lab)
% function Mesh = Mesh3D(NURBS, Lab)
% modified from igafem package (by Vinh Phu Nguyen)

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

if strcmp(Lab, 'ScalarField')
    Dof = 1; % number of dof per control point
elseif strcmp(Lab, 'VectorField')
    Dof = 3;
else
    error('"lab" must be "ScalarField" or "VectorField"')
end

NElDir = zeros(size(NURBS.KntVect)); % number of elements per direction
for i = 1 : NURBS.Dim
    NElDir(i) = numel(NURBS.uqKntVect{i}) - 1;
end

% chan =
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
%        13 14 15 16
% for a 4x2x2 control points


chan = zeros(NURBS.NCtrlPts);
chan = permute(chan, [3 2 1]);

count = 1;

for i = 1 : NURBS.NCtrlPts(3)
    for j = 1 : NURBS.NCtrlPts(2)
        for k = 1 : NURBS.NCtrlPts(1)
            chan(i, j, k) = count;
            count = count + 1;
        end
    end
end

% determine our element ranges and the corresponding
% knot Indices along each direction

ElConn{1} = BuildConn(NURBS.Order(1), NURBS.KntVect{1});
ElConn{2} = BuildConn(NURBS.Order(2), NURBS.KntVect{2});
ElConn{3} = BuildConn(NURBS.Order(3), NURBS.KntVect{3});

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4
% for a 4x2 mesh

NEl = prod(NElDir); % number of elements
NEN = prod(NURBS.Order + 1); % number of local basis functions
El = zeros(NEl, NEN);

e = 1;
for iZeta = 1 : NElDir(3)
    conn{3} = ElConn{3}(iZeta, :);
    for iEta = 1 : NElDir(2)
        conn{2} = ElConn{2}(iEta, :);
        for ixi = 1 : NElDir(1)
            c = 1;
            conn{1} = ElConn{1}(ixi, :);
            for i = 1 : numel(conn{3})
                for j = 1 : numel(conn{2})
                    for k = 1 : numel(conn{1})
                        El(e, c) = chan(conn{3}(i), conn{2}(j), conn{1}(k));
                        c = c + 1;
                    end
                end
            end
            e = e + 1;
        end
    end
end

Mesh.Dof = Dof;
Mesh.El = El;
Mesh.NEN = NEN;
Mesh.NEl = NEl;
Mesh.NElDir = NElDir;
Mesh.ElConn = ElConn;
Mesh.NDof = NURBS.NNP * Dof; % scalar field

for iSide = 1 : 6
    ind = setdiff(1 : 3, ceil(iSide / 2)); %ind=[2 3; 2 3; 1 3; 1 3; 1 2; 1 2]
    ind2 = floor((iSide + 1) / 2); % ind2 = [1 1 2 2 3 3];
    
    NCtrlPtsDir = NURBS.NCtrlPts(ind); % number of control points per direction
    Boundary.NDofDir = NCtrlPtsDir;
    NCtrlPtsBdry = prod(NCtrlPtsDir);
    Boundary.NDof = NCtrlPtsBdry * Dof;
    
    idx = ones (3, NCtrlPtsBdry);
    [idx(ind(1),:), idx(ind(2),:)] = ind2sub(NCtrlPtsDir, 1 : NCtrlPtsBdry);
    if (rem (iSide, 2) == 0)
        idx(ind2,:) = NURBS.NCtrlPts(ind2);
    end
    Boundary.Dofs = sub2ind(NURBS.NCtrlPts, idx(1, :), idx(2, :), idx(3, :));
    
    if Dof == 3 % vector field
        UDofs = sub2ind (NURBS.NCtrlPts, idx(1,:), idx(2,:), idx(3,:));
        
        Boundary.Dofs = [UDofs, UDofs + NURBS.NNP, UDofs + NURBS.NNP * 2];
        
        Boundary.CompDofs{1} = UDofs;
        Boundary.CompDofs{2} = UDofs + NURBS.NNP;
        Boundary.CompDofs{3} = UDofs + NURBS.NNP * 2;
    elseif Dof == 1 % scalar field
        Boundary.Dofs = sub2ind(NURBS.NCtrlPts, idx(1, :), idx(2, :), idx(3, :));
    end
    
    % connectivity for mesh Boundary
    UConn = reshape(ElConn{ind(1)}', NURBS.Order(ind(1)) + 1, 1, NElDir(ind(1)), 1);
    UConn = repmat(UConn, [1, NURBS.Order(ind(2)) + 1, 1, NElDir(ind(2))]);
    UConn = reshape(UConn, [], prod(NElDir(ind)));
    
    VConn = reshape(ElConn{ind(2)}', 1, NURBS.Order(ind(2)) + 1, 1, NElDir(ind(2)));
    VConn = repmat(VConn, [NURBS.Order(ind(1)) + 1, 1, NElDir(ind(1)), 1]);
    VConn = reshape(VConn, [], prod(NElDir(ind)));
    
    ElConnBdry = zeros(prod(NURBS.Order(ind) + 1), prod(NElDir(ind)));
    Indices = (UConn ~= 0) & (VConn ~= 0);
    ElConnBdry(Indices) = ...
        sub2ind(NURBS.NCtrlPts(ind), UConn(Indices), VConn(Indices));
    Boundary.El = ElConnBdry';
    Boundary.NEl = prod(NElDir(ind));
    Boundary.NElDir = NElDir(ind);
    Boundary.NEN =  prod(NURBS.Order(ind) + 1);
    Mesh.Boundary(iSide) = Boundary;
    clear Boundary
end
end




