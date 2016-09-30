function [R0, R1, R2] = RationalizeBeam(W, N0, N1, N2)
% function [R0, R1] = Rationalize(W, N0, N1)
% ----------------------------------------------------------------------
% Rationalize Bspline basis functions and their first derivatives
% ----------------------------------------------------------------------
% Input:
%       W: weights of control points, size(W) = [1, NEN]
%       N0: bspline basis functions, size(N0) = [1, NEN]
%       N1: first derivatives of bsplines, size(N1) = [Dim, NEN],
%           Dim is dimension of problem (Dim \neq NSD)
%-----------------------------------------------------------------------
% Output:
%       R0: NURBS basis functions
%       R1: corresponding first derivatives
% ----------------------------------------------------------------------

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

R0 = N0 .* W;
W0 = sum(R0);
R0 = R0 / W0;
R1 = zeros(size(N1));
R2 = zeros(size(N2));
for i = 1 : size(N1, 1)
    R1(i, :) = N1(i, :) .* W - R0 * sum(N1(i, :) .* W);
end
R1 = R1 / W0;
for i = 1 : size(N2, 1)
    R2(i, :) = N2(i, :) .* W - 2*sum(N1(i, :) .* W)*R1(i, :) - R0 * sum(N2(i, :) .* W);
end
R2 = R2 / W0;
end