function [f, d, FreeIdcs] = applyDrchltBdryVals(BdryIdcs, BdryVals, K, f)
% function [f, d, freeIdcs] = applyDrchltBdryVals(BdryIdcs, BdryVals, K, f)
% Apply boundary values for Dirichlet boundary condition
% ----------------------------------------------------------------
% Input: 
%       BdryIdcs: indices of boundary control points
%       BdryVals: corresponding values at these control points
%       K: global stiffness matrix
%       f: force vector (right hand side)
% Output:
%       f: modified force vector
%       d: displacement vector
%       freeIdcs: vector contains indices of unknowns
% -----------------------------------------------------------------

NDof = numel(f);
FreeIdcs = setdiff(1 : NDof, BdryIdcs);
d = zeros(NDof, 1);
d(BdryIdcs) = BdryVals;
f(FreeIdcs) = f(FreeIdcs) - K(FreeIdcs, BdryIdcs) * BdryVals; 
end