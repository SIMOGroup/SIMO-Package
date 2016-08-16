function [J2, W, Xi, N] = calcDersBasisFunsAtGPs(p, n, KntVect, d, NGPs, NEl)
% [J2, W, Xi, N] = calcDersBasisFunsAtGPs(p, n, KntVect, d, NGPs, NEl)
% Evaluate univariate basis functions and corresponding derivatives
% at Gauss points
%-------------------------------------------------------------
% Input:
%       p: degree of basis functions
%       n: number of control points
%       KntVect: knot vector
%       d: degree of derivative need to compute
%       NGPs: number of gauss points
%       NEl: number of elements
% ------------------------------------------------------------
% Output:
%       J2: jacobian of mapping from parent element to parametric space
%       W: gaussian quadrature weights
%       X: coordinate of Gauss points in parametric space
%       N: basis functions and their derivatives
% -------------------------------------------------------------------

Idx = zeros(1, NEl); % span index
iE = 1;
for i = p : n
    if (abs(KntVect(i) - KntVect(i + 1)) > eps(KntVect))
        Idx(iE) = i;
        iE = iE + 1;
    end
end
[Xg, W] = GaussRule(NGPs);

J2 = zeros(1, NEl);
Xi = zeros(NEl, NGPs);
N = zeros(NEl, NGPs, p + 1, d + 1);
for ex = 1 : NEl
    i = Idx(ex);
    % jacobian of mapping from parent element to parametric space
    J2(ex) = (KntVect(i + 1) - KntVect(i)) / 2;
    % compute coordinate in parameter space
    Xi(ex, :) = (Xg + 1) * J2(ex) + KntVect(i);
    % compute element basis functions, first derivative of
    % basis functions w.r.t parametric coordinate
    N(ex, :, :, :) = DersBasisFuns(repmat(i, 1, NGPs), Xi(ex, :), p,...
        d, KntVect);    
end
end