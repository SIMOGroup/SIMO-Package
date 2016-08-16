function idx = FindSpan(n, p, pts, KntVect)
% idx = FindSpan(n, p, pts, Knts)
%----------------------------------------------------------
% Determine the knot span index
% ie. find i if knot(j) lies in the half-open interval 
% [knot(i), knot(i + 1)), i = 1, 2,...,n.
%----------------------------------------------------------
% Input:
%    n: number of control points (or basis functions)
%    p: order of Bspline basis functions
%    pts: parametric points
%    KntVect: knot vector
%----------------------------------------------------------
% Output:
%    idx: knot span index
%----------------------------------------------------------
% Based on Algorithm A2.1 [The NURBS BOOK, p.68]
%----------------------------------------------------------

idx = zeros(size(pts));
for i = 1 : numel(pts)
    if (pts(i) == KntVect(n + 1))
        idx(i) = n;
        return
    end
    low = p + 1; %KntVect={a,...,a,knot_{p+2},...,knot_{n},b,...,b}
    high = n + 1;
    mid = floor((low + high)/2);
    while(pts(i) < KntVect(mid) ||...
            pts(i) >= KntVect(mid + 1))
        if(pts(i) < KntVect(mid))
            high = mid;
        else
            low = mid;
        end
        mid = floor((low + high) / 2);
    end
    idx(i) = mid;
end
end