function El = BuildConn(p, KntVect)
% El = BuildConn(p, KntVect)
% compute connectivity matrix for 1D NURBS
% -------------------------------------------------------------------------
% Input:
%       p: degree of basis function
%       KntVect: knot vector
% -------------------------------------------------------------------------
% Output:
%       El: connectivity matrix, size(El) = [NEL x NEN]
% -------------------------------------------------------------------------
% modified from igafem code of Vinh Phu Nguyen

NEl = numel(unique(KntVect)) - 1;
El = zeros(NEl, p + 1);

KntIdx = zeros(NEl, 2);
iE = 1;
for i = p : numel(KntVect) - p - 1
    if (abs(KntVect(i) - KntVect(i + 1)) > eps(KntVect))
        KntIdx(iE, :) = [i i + 1];
        iE = iE + 1;
    end
end

NRepKnts = 0; % number of repeated knots

for iE = 1 : NEl
    indices = (KntIdx(iE, 1) - p + 1) : KntIdx(iE, 1);
    preKnts = KntVect(indices);
    curKnts = ones(1, p) * KntVect(KntIdx(iE, 1));
    if isequal(preKnts, curKnts) && length(nonzeros(preKnts)) > 1;
        NRepKnts = NRepKnts + 1;
    end
    El(iE, :) = (KntIdx(iE, 1) - p) : KntIdx(iE, 1);
end
end
