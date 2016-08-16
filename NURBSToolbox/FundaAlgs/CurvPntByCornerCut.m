function Cw = CurvPntByCornerCut(n, p, KntVect, Pw, Knt)
% Cw = CurvPntByCornerCut(n, p, KntVect, Pw, Knt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute point on rational B-spline curve
% ---------------------------------------------------------------
% Input:
%       n: number of basis functions
%       p: order (degree) of basis functions
%       KntVect: knot vector
%       Pw: control points
%       Knt: a knot value
% ---------------------------------------------------------------
% Output:
%       Cw: new control points
% ---------------------------------------------------------------
% Based on Algorithm A5.2 [The NURBS BOOK, p.153].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Knt == KntVect(1)
    Cw = Pw(:, 1);
elseif Knt == KntVect(n + p + 1)
    Cw = Pw(:, n);
else    
    k = FindSpan(n, p, Knt, KntVect);
    s = FindMult(Knt, KntVect);
    r = p - s;
    Rw(:, 1 : r + 1) = Pw(:, k - p : k - p + r);
    for j = 1 : r
        for i = 1 : r - j + 1
            alfa = (Knt - KntVect(k - p + j + i - 1)) /...
                (KntVect(i + k) - KntVect(k - p + j + i - 1));
            Rw(:, i) = alfa * Rw(:, i + 1) + (1 - alfa) * Rw(:, i);
        end
    end
    Cw = Rw(:, 1);
end
end