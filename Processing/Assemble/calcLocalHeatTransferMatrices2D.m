function [HVals, FVals] = calcLocalHeatTransferMatrices2D(NURBS, Mesh, h, Ta, Ref)
% function [HVals, FVals] = calcLocalHeatTransferMatrices2D(NURBS, Mesh, h, Ta, Ref)
% -----------------------------------------------------------------
% Evaluate local surface heat transfer matrices for 2D
% -----------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       Mesh: Mesh structure
%       h: convective heat transfer coefficient
%       Ta: ambient temperature
%       Ref: referenced index to indicate boundary curve,
%       ******************************************************
%       *                                                    *
%       *                   Refs(i) = 4                      *
%       *                   (v = 1)                          *
%       *               o------------o              v        *
%       *               |            |              ^        *
%       *               |            |              |        *
%       *   Refs(i) = 1 |            | Refs(i) = 2  |        *
%       *       (u = 0) |            |   (u = 1)    +----> u *
%       *               |            |                       *
%       *               o------------o                       *
%       *                  Refs(i) = 3                       *
%       *                   (v = 0)                          *
%       ******************************************************
% --------------------------------------------------------------------
% Output:
%       HVals: matrix stores surface heat transfer matrices
%           (size(HVals) = [NEN ^ 2, NEL])
%       FVals: matrix stores local body forces
%           (size(FVals) = [NEN, NEL])
% --------------------------------------------------------------------

MeshBdry = Mesh.Boundary(Ref);
NURBSBdry = NURBSBoundary(NURBS, Ref);
NGPs = NURBSBdry.Order + 1;
[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBSBdry.Order, NURBSBdry.NCtrlPts, NURBSBdry.KntVect{1}, 1, NGPs, MeshBdry.NElDir);
CtrlPts = NURBSBdry.CtrlPts3D(1 : 2, :)';
Weights = NURBSBdry.Weights;

HVals = zeros(MeshBdry.NEN ^ 2, MeshBdry.NEl);
FVals = zeros(MeshBdry.NEN, MeshBdry.NEl);
for ex = 1 : MeshBdry.NElDir
    HL = zeros(MeshBdry.NEN); % local surface heat transfer matrix
    FL = zeros(MeshBdry.NEN, 1);
    for qx = 1 : NGPs
        
        N0 = reshape(Nx(ex, qx, :, 1), 1, []);
        N1 = reshape(Nx(ex, qx, :, 2), 1, []);
        
        [R0, R1] = Rationalize(Weights(MeshBdry.El(ex, :)), N0, N1);
        
        dxdxi = R1 * CtrlPts(MeshBdry.El(ex, :), :);
        
        J1 = norm(dxdxi);
        
        HL = HL + R0' * h * R0 * J1 * Jx(ex) * Wx(qx);
        
        FL = FL + R0' * h * Ta * J1 * Jx(ex) * Wx(qx);
    end
    FVals(:, ex) = FL;
    HVals(:, ex) = HL(:);
end
end