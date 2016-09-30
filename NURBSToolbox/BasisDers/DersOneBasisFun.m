function N0n = DersOneBasisFun(p, KntVect, ith, Pts, n)
%--------------------------------------------------------------------------
% Compute derivatives of basis function Nip
%-------------------------------------------------------------
% Input:
%       p: order of basis function
%       KntVect: knot vector
%       ith: the ith basis function
%       Pts: parametric points
%       n: maximum order of derivatives
% ------------------------------------------------------------
% Output:
%       N0n: basis functions and their derivatives
% ------------------------------------------------------------
% Based on Algorithm A2.5 [The NURBS BOOK, p.76].
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

N0n = zeros(numel(Pts), n + 1);

for i =  1 : numel(Pts)    
    iXi = Pts(i);    
    Ni = zeros(1, n + 1);
    N = zeros(p + 1, n + 1);   
    if iXi < KntVect(ith) || iXi > KntVect(ith + p + 1)
        continue;
    end    
    if iXi == KntVect(ith + p + 1)
        iXi = KntVect(ith + p + 1) - 1e-6;
    end    
    for j = 0 : p % Initialize zeroth-degree functions
        if ((iXi >= KntVect(ith + j)) &&...
                (iXi < KntVect(ith + j + 1)))
            N(j + 1, 1) = 1;
        else
            N(j + 1, 1) = 0;
        end
    end
    % Compute full triangular table
    for k = 1 : p
        if (N(1, k) == 0)
            saved = 0;
        else
            saved = ((iXi - KntVect(ith)) * N(1, k)) /...
                (KntVect(ith + k) - KntVect(ith));
        end
        for j = 0 : p - k
            Xi_left = KntVect(ith + j + 1);
            Xi_right = KntVect(ith + j + k + 1);
            if (N(j + 2, k) == 0)
                N(j + 1, k + 1) = saved;
                saved = 0;
            else
                temp = N(j + 2, k) / (Xi_right - Xi_left);
                N(j + 1, k + 1) = saved + (Xi_right - iXi) * temp;
                saved = (iXi - Xi_left) * temp;
            end
        end
    end
    Ni(1) = N(1, p + 1); % The function value    
    ND = zeros(1, p + 1);
    for k = 1 : n % Compute the derivatives
        for j = 0 : k % Load appropriate column
            ND(j + 1) = N(j + 1, p - k + 1);
        end
        for jj = 1 : k % Compute table of width k
            if ND(1) == 0
                saved = 0;
            else
                saved = ND(1) / (KntVect(ith + p - k + jj)...
                    - KntVect(ith));
            end
            for j = 0 : k - jj
                Xi_left = KntVect(ith + j + 1);
                Xi_right = KntVect(ith + j + p - k + jj + 1);
                if ND(j + 2) == 0
                    ND(j + 1) = (p - k + jj) * saved;
                    saved = 0;
                else
                    temp = ND(j + 2) / (Xi_right - Xi_left);
                    ND(j + 1) = (p - k + jj) * (saved - temp);
                    saved = temp;
                end
            end
        end
        Ni(k + 1) = ND(1);
    end
    N0n(i, :) = Ni;
end
end