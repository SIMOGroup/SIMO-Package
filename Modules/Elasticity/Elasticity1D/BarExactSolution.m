function solu = BarExactSolution(q1, q2, E, A, L, u0, P)
% function solu = BarExactSolution(q1, q2, E, A, L)
% Calculate exact solution for bar problem
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

solu.displ = @(x) (q1 - q2) * x .^ 3 / (6 * E(x) * A(x) * L) -...
    q1 * x .^ 2 / (2 * E(x) * A(x)) + (3 * L ^ 2 * q1 + ...
    3 * L ^ 2 * q2 + 6 * L * P) * x / (6 * E(x) * A(x) * L) + u0;
solu.sigma = @(x) ((q1 - q2) * x .^ 2 / (2 * E(x) * A(x) * L) -...
    q1 * x / (E(x) * A(x)) + (3 * L ^ 2 * q1 + ...
    3 * L ^ 2 * q2 + 6 * L * P) / (6 * E(x) * A(x) * L)) * E(x);
end