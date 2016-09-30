function solu = CantiliverBeamExactSolu(L, D, E, I, nu, P)
% function solu = CantiliverBeamExactSolu(L, D, E, I, nu, P)
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

solu.ux = @(x, y) P * y / (6 * E * I) .* ((6 * L - 3 * x) .* x ... 
          + (2 + nu) * (y .^ 2 - D ^ 2 / 4));
solu.uy = @(x, y) - P / (6 * E * I) * (3 * nu * y .^ 2 .* (L - x)... 
          + (4 + 5 * nu) * D ^ 2 .* x / 4 + (3 * L - x) .* x .^ 2);
end