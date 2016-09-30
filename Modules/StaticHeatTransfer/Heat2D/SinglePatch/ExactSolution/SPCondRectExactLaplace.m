function T = SPCondRectExactLaplace(fx, x, y, a, b, n)
% T = SPCondRectExactLaplace(fx, x, y, a, b, n)
% Evaluate Laplace solution
% ---------------------------------------------------------------
% Input:
%       fx: function on boundary y = b
%       x: x value, 0 <= x <= a
%       y: y value, 0 <= y <= b
%       a: first dimension of the rectangular domain
%       b: second dimension of the rectangular domain
%       n: number of summation
% -------------------------------------------------------------
% Output:
%       T: solution of laplace equation at given values
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

nArray = 1 : n;
fun = @(x) fx(x) * sin((nArray * pi * x) / a);
bn = (2 / a * integral(fun, 0, a, 'ArrayValued',true)) ./ (sinh((nArray * pi * b) / a));
tmp = zeros([size(x), n]);
for i = 1 : n
    tmp(:, :, i) = bn(i) * sin(i * pi * x/a) .* sinh(i * pi * y/a);
end
T = sum(tmp, 3);
end