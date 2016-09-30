function [C, Stress] = StressEval(NURBS, d, ParaPts, E, nu, lab)
% [C, Stress] = StressEval(NURBS, d, ParaPts, E, nu, lab)

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

D = getElastMat(E, nu, lab);
F1 = reshape(d, NURBS.NNP, [])';
F2 = reshape(reshape(d, NURBS.NNP, [])', [size(F1, 1), NURBS.NCtrlPts]);
g = GradEval(NURBS.KntVect, NURBS.CtrlPts4D, ParaPts, F2);
NPts = cellfun(@numel, ParaPts);
Cw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, ParaPts);
if NURBS.Dim == 2
    Stress = zeros([3, NPts]);
    for i = 1 : NPts(1)
        for j = 1 : NPts(2)
            Strain(1) = g(1, 1, i, j); % \epsilon_{xx}
            Strain(2) = g(2, 2, i, j); % \epsilon_{yy}
            Strain(3) = g(2, 1, i, j) + g(1, 2, i, j); % 2\epsilon_{xy}
            Stress(:, i, j) = D * Strain';
        end
    end
    Weights = Cw(4, :, :);
    C = bsxfun(@rdivide, Cw(1 : 3, :, :), Weights); 
end
end