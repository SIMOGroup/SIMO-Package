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

%{
Copyright (C) <2014-2016>  <Khanh Chau-Nguyen, Hung Nguyen-Xuan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

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