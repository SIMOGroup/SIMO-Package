function Ni = OneBasisFun(p, KntVect, i, pts)
% Ni = OneBasisFun(p, KntVect, i, pts)
%--------------------------------------------------------------------------
% Compute the ith basis function
%---------------------------------------------------------------
% Input:
%       p: order of basis function
%       KntVect: knot vector
%       i: the ith basis function
%       pts: parametric points
%---------------------------------------------------------------
% Output:
%       Ni: individual B-spline basis function
%---------------------------------------------------------------
% Based on Algorithm A2.4 [The NURBS BOOK, p.74]
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

m = numel(KntVect);
Ni = zeros(numel(pts), 1);
for ii = 1 : numel(pts)
    xi_i = pts(ii);
    % Special case
    if(i == 1 && xi_i == KntVect(1) || i == m - p - 1 &&...
            xi_i == KntVect(m))
        Nip = 1;
        Ni(ii) = Nip;
        continue;
    end
    % Local property
    if(xi_i < KntVect(i) || xi_i >= KntVect(i + p + 1))
        Nip = 0;
        Ni(ii) = Nip;
        continue;
    end
    N = zeros(1, p + 1);
    % Initialize the zeroth - degree functions
    for j = 0 : p
        if ((xi_i >= KntVect(i+j)) &&...
                (xi_i < KntVect(i + j + 1)))
            N(j + 1) = 1;
        else
            N(j + 1) = 0;
        end
    end
    % Compute triangular table
    for k = 1 : p
        if (N(1) == 0)
            saved = 0;
        else
            saved = ((xi_i - KntVect(i)) * N(1)) /...
                (KntVect(i + k) - KntVect(i));
        end
        for j = 0 : p - k
            Xi_left = KntVect(i + j + 1);
            Xi_right = KntVect(i + j + k + 1);
            if (N(j + 2) == 0)
                N(j + 1) = saved;
                saved = 0;
            else
                temp = N(j + 2) / (Xi_right - Xi_left);
                N(j + 1) = saved + (Xi_right - xi_i) * temp;
                saved = (xi_i - Xi_left) * temp;
            end
        end
        Nip = N(1);
    end
    Ni(ii) = Nip;
end
end