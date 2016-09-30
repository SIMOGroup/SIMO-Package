function u2 = SPCondRectExactPoisson(s, xVal, yVal, aVal, bVal, nVal, mVal)
% function u2 = SPCondRectExactPoisson(s, xVal, yVal, aVal, bVal, nVal, mVal)
% Evaluate homogenous Poisson solution
% ---------------------------------------------------------------
% Input:
%       s: volumetric source function (W m^{-2})
%       xVal: x value, 0 \leq x \leq a
%       yVal: y value, 0 \leq y \leq b
%       aVal: first dimension of the rectangular domain
%       bVal: second dimension of the rectangular domain
%       nVal: number of summation in first direction
%       mVal: number of summation in second direction
% -------------------------------------------------------------
% Output:
%       u2: solution of Poisson equation at given values
% -------------------------------------------------------------

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

syms x y a b m n
% fourier coefficients
nmrtr = 4 * (int(int(s * sin((n * pi * x) / a) * sin((m * pi * y) / b), x, [0, a]), y, [0, b] ));
dnmntr = (pi ^ 2) * a * b * ((n ^ 2)/ (a ^ 2) + (m ^ 2 ) / (b ^ 2));
anm = nmrtr / dnmntr;
% fourier series terms
uxy = anm * sin((n * pi * x) / a) * sin((m * pi * y) / b);
% solution of homogenous poisson problem
% symbolic summation
uxySum = symsum(symsum(uxy, n, 1, nVal), m, 1, mVal);
% evaluate the solution at the given input values
u2 = eval(subs(uxySum, {x, y, b, a}, {xVal, yVal, bVal, aVal}));
end