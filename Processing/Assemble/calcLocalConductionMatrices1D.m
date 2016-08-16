function [KVals, FVals] = calcLocalConductionMatrices1D(Mesh, NURBS, k, s, LAB, varargin)
% [KVals, FVals] = calcLocalConductionMatrices1D(Mesh, NURBS, k, s, LAB, varargin)
% --------------------------------------------------------------------
% Evaluate local conduction matrices and body forces
% --------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       ka: conductivity coefficient
%       s: heat source
%       LAB: LABEL keyword to indicate the problem type, it can be 'FIN' 
%       for cooling fin problem or 'PLANE WALL' for normal problem
%       varargin: parameters for cooling fin problem, where
%           varargin{1} = P: perimeter of cross section
%           varargin{2} = A: area of cross section
%           varargin{3} = h: convective heat transfer coefficient.
%           varargin{4} = Ta: ambient temperature
%       i.e.
%       (Mesh, NURBS, k, s, 'FIN', P, A, h, Ta)
% --------------------------------------------------------------------
% Output:
%       KVals: matrix stores local conduction matrices
%           (size(KVals) = [NEN ^ 2, NEL])
%       FVals: matrix stores local body forces
%           (size(FVals) = [NEN, NEL])
% --------------------------------------------------------------------

NGPs = NURBS.Order + 1; % number of gauss points

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order, NURBS.NCtrlPts(1), NURBS.KntVect{1}, 1, NGPs, Mesh.NEl);

KVals = zeros(Mesh.NEN ^ 2, Mesh.NEl);
FVals = zeros(Mesh.NEN, Mesh.NEl);

Weights = NURBS.Weights;
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';
if strcmp(LAB, 'PLANE WALL')
    A = 1;
    P = 0;
    Ta = 0;
    h = 0;
elseif strcmp(LAB, 'FIN')
    P = varargin{1}; % perimeter
    A = varargin{2}; % area
    h = varargin{3};
    Ta = varargin{4};
else
    error('You must specify the LABEL keyword')
end

% Loop over elements (knot spans)
for e = 1 : Mesh.NEl
    Ke = zeros(Mesh.NEN);
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
        dRdx = J1 \ R1;
        
        % compute elementary stiffness matrix
        Ke = Ke + (dRdx' * k * A * dRdx + R0' * h * P * R0) * J1 * Jx(e) * Wx(qx);
        
        %         % compute distributed force
        Pts = R0 * CtrlPts(Mesh.El(e, :), :); % physical coords of gauss points
        
        x = Pts(1);
        
        Fe = Fe + R0' * (P * h * Ta + s(x) * A) * J1 * Jx(e) * Wx(qx);
    end
    FVals(:, e) = Fe;
    KVals(:, e) = Ke(:);
end
end