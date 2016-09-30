function solu = BeamExactSolution(L, k, q, F)
% function solu = BeamExactSolution(E, I, L, k, q, F)
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

solu.displ = @(x)F/(6*k)*(3*L*x.^2-x.^3);
%solu.moment = @(x)q/(24*k)*(x.^4 - 4*L*x.^3 + 6*L^2*x.^2); for
%distribution
end