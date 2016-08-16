function [KVals, FVals] = calcLocalStiffnessMatrices1D(Mesh, NURBS, E, A, df)
% function [KVals, FVals] = calcLocalStiffnessMatrices1D(Mesh, NURBS, E, A, df)
% ----------------------------------------------------------------------
% Calculate one dimensional local stiffness matrices
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure of the computational domain
%       NURBS: NURBS structure
%       E: Young's modulus
%       A: area of cross section
%       df: function of distributed force
% ---------------------------------------------------------------------
% Output:
%       KVals: matrix stores local stiffness matrices
%           (size(KVals) = [NEN ^ 2, NEL])
%       FVals: matrix stores local body forces
%           (size(FVals) = [NEN, NEL])
% --------------------------------------------------------------------

NGPs = NURBS.Order(1) + 1; % number of gauss points

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs, Mesh.NEl);

KVals = zeros(Mesh.NEN ^ 2, Mesh.NEl);
FVals = zeros(Mesh.NEN, Mesh.NEl);

Weights = NURBS.Weights;
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';
% Loop over elements (knot spans)
for e = 1 : Mesh.NEl
    Ke = zeros(Mesh.NEN, Mesh.NEN);
    Fe = zeros(Mesh.NEN, 1);
    % loop over Gauss points
    for qx = 1 : NGPs
        
        N0 = Nx(e, qx, :, 1); % bspline basis functions
        N1 = Nx(e, qx, :, 2); % 1st derivarive of bspline
        
        % rationalize to create nurbs basis functions and their derivatives
        [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0(:)', N1(:)');
        
        % gradient of mapping from parameter space to physical space
        dxdxi = R1 * CtrlPts(Mesh.El(e, :), :);
        % compute the jacobian of physical and parameter domain mapping
        J1 = norm(dxdxi);
        
        % compute derivative of basis functions w.r.t spatial
        % physical coordinates
        dRdx = J1 ^ (-1) * R1;
        
        % physical coords of gauss points
        Pts = R0 * CtrlPts(Mesh.El(e, :), :);         
        x = Pts(1);   
        
        % compute local stiffness matrix
        Ke = Ke + dRdx' * A(x) * E(x) * dRdx * J1 * Jx(e) * Wx(qx);
        
        % compute distributed force     
        Fe = Fe + df(x) * R0' * J1 * Jx(e) * Wx(qx);
    end
    FVals(:, e) = Fe;
    KVals(:, e) = Ke(:);
end
end