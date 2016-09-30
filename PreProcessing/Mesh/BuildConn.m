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

NEl = numel(unique(KntVect)) - 1;
El = zeros(NEl, p + 1);

KntIdx = zeros(NEl, 2);
iE = 1;
for i = p : numel(KntVect) - p - 1
    if (abs(KntVect(i) - KntVect(i + 1)) > sqrt(eps))
        KntIdx(iE, :) = [i, i + 1];
        iE = iE + 1;
    end
end

for iE = 1 : NEl
    El(iE, :) = (KntIdx(iE, 1) - p) : KntIdx(iE, 1);
end
end
