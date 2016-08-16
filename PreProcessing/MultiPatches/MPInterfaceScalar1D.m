function [GNum, GNDof] = MPInterfaceScalar1D(Interfaces, Mesh)
% function [GNum, GNDof] = MPInterfaceScalar1D(Interfaces, Mesh)

if (~isempty(Interfaces))
    GNum = cell(numel(Mesh), 1);
    patch_intrfc = cell(numel(Mesh), 1);
    ttform = cell(numel(Mesh), numel(Interfaces));
    ppnum = cell(numel(Interfaces), 1);
    for iPtc = 1:numel(Mesh)
        GNum{iPtc} = zeros (1, Mesh{iPtc}.NDof);
        patch_intrfc{iPtc} = union (find([Interfaces.Patch1] == iPtc), ...
            find([Interfaces.Patch2] == iPtc));
    end
    
    for intrfc = 1:numel(Interfaces)
        iptc1  = Interfaces(intrfc).Patch1;
        iptc2  = Interfaces(intrfc).Patch2;
        iside1 = Interfaces(intrfc).Side1;
        iside2 = Interfaces(intrfc).Side2;
        ttform{iptc1, intrfc} = Mesh{iptc1}.Boundary(iside1).Dofs;        
        ttform{iptc2, intrfc} = Mesh{iptc2}.Boundary(iside2).Dofs;
        ppnum{intrfc} = 0;
    end
    
    GNDof = 0;
    % We start with the dofs that do not belong to any interface
    for iPtc = 1:numel (Mesh)
        non_intrfc_dofs = setdiff(1:Mesh{iPtc}.NDof, [ttform{iPtc,:}]);
        GNum{iPtc}(non_intrfc_dofs) = GNDof + (1:numel(non_intrfc_dofs));
        GNDof = GNDof + numel (non_intrfc_dofs);
    end
    
    % Then we set the interfaces
    for intrfc = 1:numel(Interfaces)
        iPtc     = Interfaces(intrfc).Patch1;
        new_dofs = find (GNum{iPtc}(ttform{iPtc, intrfc}) == 0);
        new_dofs_number = GNDof + (1:numel (new_dofs));
        GNum{iPtc}(ttform{iPtc, intrfc}(new_dofs)) = new_dofs_number;
        ppnum{intrfc}(new_dofs) = new_dofs_number;
        
        [GNum, ppnum] = set_same_interface (iPtc, intrfc, ttform, new_dofs,...
            Interfaces, GNum, ppnum, patch_intrfc);
        [GNum, ppnum] = set_same_patch (iPtc, intrfc, ttform, new_dofs, ...
            Interfaces, GNum, ppnum, patch_intrfc);
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

function [glob_num, ppnum] = set_same_patch (iptc, intrfc, ttform, new_dofs, ...
    interfaces, glob_num, ppnum, patch_intrfc)

intrfc_dofs = ttform{iptc, intrfc}(new_dofs);
for ii = setdiff (patch_intrfc{iptc}, intrfc)
    [common_dofs, pos1, pos2] = intersect (ttform{iptc, ii}, intrfc_dofs);
    not_set = find (ppnum{ii}(pos1) == 0);
    if (~isempty (not_set))
        ppnum{ii}(pos1(not_set)) = ppnum{intrfc}(new_dofs(pos2(not_set)));
        [glob_num, ppnum] = set_same_interface (iptc, ii, ttform, ...
            pos1(not_set), interfaces, glob_num, ppnum, patch_intrfc);
    end
end
end

function [glob_num, ppnum] = set_same_interface (iptc, intrfc, ttform, ...
    new_dofs, interfaces, glob_num, ppnum, patch_intrfc)
intrfc_dofs = ttform{iptc, intrfc}(new_dofs);
iptc2 = setdiff ([interfaces(intrfc).Patch1 interfaces(intrfc).Patch2], iptc);
glob_num{iptc2}(ttform{iptc2, intrfc}(new_dofs)) = ...
    glob_num{iptc}(intrfc_dofs);
[glob_num, ppnum] = set_same_patch (iptc2, intrfc, ttform, new_dofs, ...
    interfaces, glob_num, ppnum, patch_intrfc);

end
