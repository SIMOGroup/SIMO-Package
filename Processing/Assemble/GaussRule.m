function [X, W] = GaussRule(n)
% function [X, W] = GaussRule(n)
% -------------------------------------------------------------
% Gauss-Legendre quadrature by Davis-Rabinowitz method.
% -------------------------------------------------------------
% The integral:
%      \int_{-1}^{1} f(x) dx
% The quadrature rule:
%      \sum_{i = 1} ^ n w(i) f(x(i)) 
%---------------------------------------------------------------
% Input:
%      n: the order
%--------------------------------------------------------------
% Output:
%      X: the abscissas
%      W: the weights
%--------------------------------------------------------------
% Reference:
%      Philip Davis, Philip Rabinowitz,
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

iter = 2; % iteration counter
e1 = n * (n + 1);
m = fix((n + 1) / 2);
t = (3 : 4 : 4 * m - 1) * pi / (4 * n + 2);
x0 = (1 - (1 - 1 / n) / (8 * n * n)) * cos(t);
for j = 1 : iter
    pkm1 = 1; 
    pk = x0;
    for k = 2 : n 
        pkp1 = 2 * x0 .* pk - pkm1 - (x0 .* pk - pkm1) / k;
        pkm1 = pk; 
        pk = pkp1;
    end
    den = 1 - x0 .* x0; 
    d1 = n * (pkm1 - x0 .* pk); 
    dpn = d1 ./ den;
    d2pn = (2 .* x0 .* dpn - e1 .* pk) ./ den;
    d3pn = (4 * x0 .* d2pn + (2 - e1) .* dpn) ./ den;
    d4pn = (6 * x0 .* d3pn + (6 - e1) .* d2pn) ./ den;
    u = pk ./ dpn; 
    v = d2pn ./ dpn;
    % Initial approximation H:
    h = -u .* (1 + 0.5 * u .* (v + u .* (v .* v - u .*...
        d3pn ./ (3*dpn))));
    % Refine H using one step of Newton's method:
    p = pk + h .* (dpn + 0.5 * h .* (d2pn + h / 3 .*...
        (d3pn + 0.25 * h .* d4pn)));
    dp = dpn + h .* (d2pn + 0.5 * h .* (d3pn + h .* d4pn / 3));
    h = h - p ./ dp;
    x0 = x0 + h;
end
X = -x0 - h;
fx = d1 - h .* e1 .* (pk + 0.5 * h .* (dpn + (h / 3) .*...
    (d2pn + 0.25 * h .* (d3pn + 0.2 * h .* d4pn))));
W = 2 * (1 - X .^ 2) ./ (fx .* fx);
if ((m + m) > n)
    X(m) = 0;
end
if (~((m + m) == n))
    m = m - 1;
end
i = 1 : m;
X(n + 1 - i) = -X(i); 
W(n + 1 - i) = W(i);
end