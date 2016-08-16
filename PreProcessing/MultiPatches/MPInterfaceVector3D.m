function [GNum, GNDof] = MPInterfaceVector3D (Interfaces, Mesh, gluedFaces)
% function [GNum, GNDof] = MPInterfaceVector3D (Interfaces, Mesh, gluedFaces)
% modified from GeoPDEs code

if (~isempty(Interfaces))
    GNum = cell(numel(Mesh), 1);
    
    patchIntrfc = cell(numel(Mesh), 1);
    ttform = cell(numel(Mesh), numel(Interfaces));
    ppnum = cell(numel(Interfaces), 1);
    for iPatch = 1 : numel(Mesh)
        GNum{iPatch} = zeros(1, Mesh{iPatch}.NDof);
        patchIntrfc{iPatch} = union(find([Interfaces.Patch1] == iPatch), find([Interfaces.Patch2] == iPatch));
    end
    
    for iInterface = 1 : numel(Interfaces)
        iPatch1  = Interfaces(iInterface).Patch1;
        iPatch2  = Interfaces(iInterface).Patch2;
        iSide1 = Interfaces(iInterface).Side1;
        iSide2 = Interfaces(iInterface).Side2;
        MeshBdry1 = Mesh{iPatch1}.Boundary(iSide1);
        ttform{iPatch1, iInterface} = MeshBdry1.Dofs;
        MeshBdry2 = Mesh{iPatch2}.Boundary(iSide2);
        nghbrDofsComp1 = reshape(MeshBdry2.CompDofs{1}, MeshBdry2.NDofDir);
        nghbrDofsComp2 = reshape(MeshBdry2.CompDofs{2}, MeshBdry2.NDofDir);
        nghbrDofsComp3 = reshape(MeshBdry2.CompDofs{3}, MeshBdry2.NDofDir);
        if (Interfaces(iInterface).Flag == -1)
            nghbrDofsComp1 = nghbrDofsComp1';
            nghbrDofsComp2 = nghbrDofsComp2';
            nghbrDofsComp3 = nghbrDofsComp3';
        end
        if (Interfaces(iInterface).Ornt1 == -1)
            nghbrDofsComp1 = flipud(nghbrDofsComp1);
            nghbrDofsComp2 = flipud(nghbrDofsComp2);
            nghbrDofsComp3 = flipud(nghbrDofsComp3);
        end
        if (Interfaces(iInterface).Ornt2 == -1)
            nghbrDofsComp1 = fliplr(nghbrDofsComp1);
            nghbrDofsComp2 = fliplr(nghbrDofsComp2);
            nghbrDofsComp3 = fliplr(nghbrDofsComp3);
        end
        nghbrDofs = [nghbrDofsComp1(:)' nghbrDofsComp2(:)' nghbrDofsComp3(:)'];
        ttform{iPatch2, iInterface} = nghbrDofs;
        ppnum{iInterface} = zeros(1, MeshBdry1.NDof);
    end
    
    GNDof = 0;
    % We start with the Dofs that do not belong to any interface
    for iPatch = 1 : numel(Mesh)
        %nonIntrfcDofs = ...
        %         setdiff(1 : Mesh{iPatch}.NDof, [ttform{iPatch,:}]);
        GluedMesh1 = Mesh{iPatch}.Boundary(gluedFaces{iPatch}(1));
        GluedMesh2 = Mesh{iPatch}.Boundary(gluedFaces{iPatch}(2));
        
        gluedDofs = union(GluedMesh1.Dofs, GluedMesh2.Dofs);
        
        nonIntrfcDofs = setdiff(1 : Mesh{iPatch}.NDof, union([ttform{iPatch, :}], gluedDofs));
        GNum{iPatch}(nonIntrfcDofs) = GNDof + (1 : numel(nonIntrfcDofs));
        GNDof = GNDof + numel(nonIntrfcDofs);
        
        dofsBdry1 = setdiff(GluedMesh1.Dofs, [ttform{iPatch, :}]);
        dofsBdry2 = setdiff(GluedMesh2.Dofs, [ttform{iPatch, :}]);
        
        newDofsBdry1 = find(GNum{iPatch}(dofsBdry1) == 0);
        newDofsBdry2 = GNum{iPatch}(dofsBdry2) == 0;
        
        newGluedDofs = GNDof + (1 : numel(newDofsBdry1));
        GNum{iPatch}(dofsBdry1(newDofsBdry1)) = newGluedDofs;
        GNum{iPatch}(dofsBdry2(newDofsBdry2)) = newGluedDofs;
        
        GNDof = GNDof + numel(newDofsBdry1);
    end
    
    % Then we set the Interfaces
    for iInterface = 1 : numel(Interfaces)
        iPatch = Interfaces(iInterface).Patch1;
        MeshBdryIntrfc1 = Mesh{iPatch}.Boundary(gluedFaces{iPatch}(1));
        MeshBdryIntrfc2 = Mesh{iPatch}.Boundary(gluedFaces{iPatch}(2));
        
        gluedDofs = union(MeshBdryIntrfc1.Dofs, MeshBdryIntrfc2.Dofs);
        
        newNonGluedDofsIntrfc = setdiff(ttform{iPatch, iInterface}, gluedDofs);
        newDofs = find(GNum{iPatch}(newNonGluedDofsIntrfc) == 0);
        newDofsNumbering = GNDof + (1 : numel (newDofs));
        GNum{iPatch}(newNonGluedDofsIntrfc(newDofs)) = newDofsNumbering;
        GNDof = GNDof + numel (newDofs);
        
        DofsBdry1Intrfc = intersect([ttform{iPatch, iInterface}], MeshBdryIntrfc1.Dofs);
        DofsBdry2Intrfc = intersect([ttform{iPatch, iInterface}], MeshBdryIntrfc2.Dofs);
        newDofsBdry1Intrfc = find(GNum{iPatch}(DofsBdry1Intrfc) == 0);
        newDofsBdry2Intrfc = GNum{iPatch}(DofsBdry2Intrfc) == 0;
        
        newGluedDofsIntrfc = GNDof + (1 : numel(newDofsBdry1Intrfc));
        GNum{iPatch}(DofsBdry1Intrfc(newDofsBdry1Intrfc)) = newGluedDofsIntrfc;
        GNum{iPatch}(DofsBdry2Intrfc(newDofsBdry2Intrfc)) = newGluedDofsIntrfc;
        
        GNDof = GNDof + numel(newDofsBdry1Intrfc);
        
        newDofs = 1 : numel(newDofs) + 2 * numel(newDofsBdry1Intrfc);
        ppnum{iInterface}(newDofs) = [newDofsNumbering, newGluedDofsIntrfc, newGluedDofsIntrfc];
        
        [GNum, ppnum] = SetSameIntrfc(iPatch, iInterface, ttform, newDofs, Interfaces, GNum, ppnum, patchIntrfc);
        [GNum, ppnum] = SetSamePatch(iPatch, iInterface, ttform, newDofs, Interfaces, GNum, ppnum, patchIntrfc);
    end
    
else
    GNDof = 0;
    GNum = cell(numel(Mesh), 1);
    for iPatch = 1 : numel(Mesh)
        GNum{iPatch} = GNDof + (1 : Mesh{iPatch}.NDof);
        GNDof = GNDof + Mesh{iPatch}.NDof;
    end
end

end

function [glob_num, ppnum] = SetSamePatch (iptc, intrfc, ttform, newDofs, Interfaces, glob_num, ppnum, patchIntrfc)

intrfcDofs = ttform{iptc, intrfc}(newDofs);
for ii = setdiff(patchIntrfc{iptc}, intrfc)
    [~, pos1, pos2] = intersect(ttform{iptc, ii}, intrfcDofs);
    notSet = find (ppnum{ii}(pos1) == 0);
    if (~isempty (notSet))
        ppnum{ii}(pos1(notSet)) = ppnum{intrfc}(pos2(notSet));
        [glob_num, ppnum] = SetSameIntrfc(iptc, ii, ttform, pos1(notSet), Interfaces, glob_num, ppnum, patchIntrfc);
    end
end
end

function [glob_num, ppnum] = SetSameIntrfc (iptc, intrfc, ttform, ...
    newDofs, Interfaces, glob_num, ppnum, patchIntrfc)
intrfcDofs = ttform{iptc, intrfc}(newDofs);
iptc2 = setdiff([Interfaces(intrfc).Patch1 Interfaces(intrfc).Patch2],iptc);
glob_num{iptc2}(ttform{iptc2, intrfc}(newDofs)) = glob_num{iptc}(intrfcDofs);
[glob_num, ppnum] = SetSamePatch (iptc2, intrfc, ttform, newDofs, Interfaces, glob_num, ppnum, patchIntrfc);
end
