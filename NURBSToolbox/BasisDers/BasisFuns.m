function N0 = BasisFuns(Idx, Pts, p, KntVect)
% N0 = BasisFuns(Idx, Pts, p, KntVect)
%--------------------------------------------------------------------------
% Compute the nonvanishing B-splines basis functions.
%---------------------------------------------------------------
% Input:
%      idx: knot span index
%      pts: parametric points
%      p: order of basis
%      KntVect: knot vector
%--------------------------------------------------------------
% Output:
%      N0: B-spline basis functions
%--------------------------------------------------------------
% Based on Algorithm A2.2 [The NURBS BOOK, p.70]
%--------------------------------------------------------------------------

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

N0 = zeros(numel(Pts), p + 1);
for i = 1 : numel(Pts)
    
    index = Idx(i);
    xi_i = Pts(i);
    
    left = zeros(1, p + 1);
    right = zeros(1, p + 1);
    
    N_i = zeros(1, p + 1);
    N_i(1) = 1;
    for j = 1 : p
        left(j + 1) = xi_i - KntVect(index + 1 - j);
        right(j + 1) = KntVect(index + j) - xi_i;
        saved = 0;
        for r = 0 : j - 1
            temp = N_i(r + 1)/(right(r + 2) + left(j - r + 1));
            N_i(r + 1) = saved + right(r + 2) * temp;
            saved = left(j - r + 1) * temp;
        end
        N_i(j + 1) = saved;
    end
    N0(i, :) = N_i;
end
end