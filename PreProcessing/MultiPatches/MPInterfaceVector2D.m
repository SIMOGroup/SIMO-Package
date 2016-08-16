% MP_INTERFACE_VECTOR_2D: create a global numbering of vectorial basis functions in two-dimensional multipatch geometries.
%
%   [GNum, GNDof] = mp_interface_vector_2d (Interfaces, sp);
%
% INPUT:
%
%  Interfaces: structure with the information of the Interfaces between patches (see mp_geo_read_file)
%  sp:         object representing the space of discrete functions (see sp_vector_2d)
%
% OUTPUT:
%
%  GNum:  global numbering of the discrete basis functions
%  GNDof: total number of degrees of freedom
%
% Copyright (C) 2011 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [GNum, GNDof] = MPInterfaceVector2D(Interfaces, Mesh)

if (~isempty (Interfaces))
    GNum = cell (numel (Mesh), 1);
    PtcIntrfc = cell(numel (Mesh), 1);
    ttform = cell(numel (Mesh), numel (Interfaces));
    ppnum = cell(numel (Interfaces), 1);
    for iPtc = 1:numel(Mesh)
        GNum{iPtc} = zeros(Mesh{iPtc}.NDof, 1);
        PtcIntrfc{iPtc} = union (find([Interfaces.Patch1] == iPtc), ...
            find([Interfaces.Patch2] == iPtc));
    end
    
    for intrfc = 1:numel(Interfaces)
        iPtc1  = Interfaces(intrfc).Patch1;
        iPtc2  = Interfaces(intrfc).Patch2;
        iSide1 = Interfaces(intrfc).Side1;
        iSide2 = Interfaces(intrfc).Side2;
        ttform{iPtc1, intrfc} = Mesh{iPtc1}.Boundary(iSide1).Dofs;
        if (Interfaces(intrfc).Ornt == 1)
            ttform{iPtc2, intrfc} = Mesh{iPtc2}.Boundary(iSide2).Dofs;
        else
            nghbrDofsComp1 = flipud(Mesh{iPtc2}.Boundary(iSide2).CompDofs{1});
            nghbrDofsComp2 = flipud(Mesh{iPtc2}.Boundary(iSide2).CompDofs{2});
            ttform{iPtc2, intrfc} = [nghbrDofsComp1; nghbrDofsComp2];
        end
        ppnum{intrfc} = zeros (1, Mesh{iPtc1}.Boundary(iSide1).NDof);
    end
    
    GNDof = 0;
    % We start with the dofs that do not belong to any interface
    for iPtc = 1 : numel(Mesh)
        non_intrfc_dofs = setdiff(1 : Mesh{iPtc}.NDof, [ttform{iPtc,:}]);
        GNum{iPtc}(non_intrfc_dofs) = GNDof + (1:numel(non_intrfc_dofs));
        GNDof = GNDof + numel (non_intrfc_dofs);
    end
    
    % Then we set the Interfaces
    for intrfc = 1:numel(Interfaces)
        iPtc     = Interfaces(intrfc).Patch1;
        new_dofs = find (GNum{iPtc}(ttform{iPtc, intrfc}) == 0);
        new_dofs_number = GNDof + (1:numel (new_dofs));
        GNum{iPtc}(ttform{iPtc, intrfc}(new_dofs)) = new_dofs_number;
        ppnum{intrfc}(new_dofs) = new_dofs_number;
        
        [GNum, ppnum] = set_same_interface (iPtc, intrfc, ttform, new_dofs,...
            Interfaces, GNum, ppnum, PtcIntrfc);
        [GNum, ppnum] = set_same_patch (iPtc, intrfc, ttform, new_dofs, ...
            Interfaces, GNum, ppnum, PtcIntrfc);
        GNDof = GNDof + numel (new_dofs);
        
    end
    
else
    GNDof = 0;
    GNum = cell (numel (Mesh), 1);
    for iPtc = 1:numel(Mesh)
        GNum{iPtc} = GNDof + (1:Mesh{iPtc}.ndof);
        GNDof = GNDof + Mesh{iPtc}.ndof;
    end
end

end

function [GNum, ppnum] = set_same_patch (iPtc, intrfc, ttform, new_dofs, ...
    Interfaces, GNum, ppnum, PtcIntrfc);

intrfc_dofs = ttform{iPtc, intrfc}(new_dofs);
for ii = setdiff (PtcIntrfc{iPtc}, intrfc)
    [common_dofs, pos1, pos2] = intersect (ttform{iPtc, ii}, intrfc_dofs);
    not_set = find (ppnum{ii}(pos1) == 0);
    if (~isempty (not_set))
        ppnum{ii}(pos1(not_set)) = ppnum{intrfc}(pos2(not_set));
        [GNum, ppnum] = set_same_interface (iPtc, ii, ttform, ...
            pos1(not_set), Interfaces, GNum, ppnum, PtcIntrfc);
    end
end
end

function [GNum, ppnum] = set_same_interface (iPtc, intrfc, ttform, ...
    new_dofs, Interfaces, GNum, ppnum, PtcIntrfc)
intrfc_dofs = ttform{iPtc, intrfc}(new_dofs);
iPtc2 = setdiff ([Interfaces(intrfc).Patch1 Interfaces(intrfc).Patch2], iPtc);
GNum{iPtc2}(ttform{iPtc2, intrfc}(new_dofs)) = ...
    GNum{iPtc}(intrfc_dofs);
[GNum, ppnum] = set_same_patch (iPtc2, intrfc, ttform, new_dofs, ...
    Interfaces, GNum, ppnum, PtcIntrfc);

end
