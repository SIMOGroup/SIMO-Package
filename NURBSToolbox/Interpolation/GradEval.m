function g = GradEval(KntVects, CtrlPts, ParaPts, d)
% g = GradEval(KntVects, CtrlPts, ParaPts, d)
% ---------------------------------------------------------------
% Evaluate gradient
% ---------------------------------------------------------------
% Input:
%       Knts: knots vector in cell format
%       CtrlPts: control points in homogenous space
%       ParaPts: parameter points to interpolate
%       d: temperature or displacement field
% --------------------------------------------------------------
% Output:
%       g: elvaluated values in matrix format
% --------------------------------------------------------------

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

assert(iscell(KntVects), 'Knot vector(s) must be in cell format');
assert(iscell(ParaPts), 'Parameter points must be in cell format');

dim = numel(KntVects);
NCtrlPts = size(CtrlPts);

p = zeros(1, dim);
Idx = cell(1, dim);
N01 = cell(1, dim);
for i = 1 : dim
    p(i) = numel(KntVects{i}) - NCtrlPts(i + 1) - 1;
    Idx{i} = FindSpan(NCtrlPts(i + 1), p(i), ParaPts{i}, KntVects{i});
    N01{i} = DersBasisFuns(Idx{i}, ParaPts{i}, p(i), 1, KntVects{i});
end

NPts = cellfun(@numel, ParaPts);

nen = prod(p + 1);
N0 = zeros(1, nen);
N1 = zeros(dim, nen);

if dim == 1
    g = zeros(size(d, 1), NPts(1));
    for i = 1 : NPts(1)
        de = d(:, Idx{1}(i) - p : Idx{1}(i));
        Weights = CtrlPts(4, Idx{1}(i) - p : Idx{1}(i));
        xe = CtrlPts(1, Idx{1}(i) - p : Idx{1}(i)) ./ Weights;
        N0 = N01{1}(i, :, 1);
        N1 = N01{1}(i, :, 2);
        [~, R1] = Rationalize(Weights, N0, N1);
        ge = R1 * de';
        dxdxi = R1 * xe';
        ge = dxdxi ^ (-1) * ge;
        g(:, i) = ge;
    end
elseif dim == 2 % surf
    g = zeros([dim, size(d, 1), NPts]);
    for i = 1 : NPts(1)
        for j = 1 : NPts(2)
            de = d(:, Idx{1}(i) - p(1) : Idx{1}(i), Idx{2}(j) - p(2) : Idx{2}(j));
            Weights = CtrlPts(4, Idx{1}(i) - p(1) : Idx{1}(i), Idx{2}(j) - p(2) : Idx{2}(j));
            xyze = bsxfun(@rdivide, CtrlPts(1 : 3, Idx{1}(i) - p(1) : Idx{1}(i), Idx{2}(j) - p(2) : Idx{2}(j)), Weights);
            
            k = 1;
            for jk = 1 : p(2) + 1
                for ik = 1 : p(1) + 1
                    N0(k) = N01{1}(i, ik, 1) * N01{2}(j, jk, 1);
                    N1(1, k) = N01{1}(i, ik, 2) * N01{2}(j, jk, 1);
                    N1(2, k) = N01{1}(i, ik, 1) * N01{2}(j, jk, 2);
                    k = k + 1;
                end
            end
            [~, R1] = Rationalize(reshape(Weights, 1, []), N0, N1);
            
            ge = R1 * reshape(de, 2, [])';
            dxdxi = R1 * reshape(xyze, [], nen)';
            
            ge = dxdxi(:, 1 : 2) ^ (-1) * ge;
            g(:, :, i, j) = ge;
        end
    end
else % volu
    g = zeros([dim, size(d, 1), NPts]);
    for i = 1 : NPts(1)
        for j = 1 : NPts(2)
            for k = 1 : NPts(3)
                de = d(:, Idx{1}(i) - p(1) : Idx{1}(i),...
                    Idx{2}(j) - p(2) : Idx{2}(j),...
                    Idx{3}(k) - p(3) : Idx{3}(k));
                
                Weights = CtrlPts(4, Idx{1}(i) - p(1) : Idx{1}(i),...
                    Idx{2}(j) - p(2) : Idx{2}(j),...
                    Idx{3}(k) - p(3) : Idx{3}(k));
                
                xyze = bsxfun(@rdivide, CtrlPts(1 : 3,...
                    Idx{1}(i) - p(1) : Idx{1}(i),...
                    Idx{2}(j) - p(2) : Idx{2}(j),...
                    Idx{3}(k) - p(3) : Idx{3}(k)), Weights);
                
                l = 1;
                for kk = 1 : p(3) + 1
                    for jk = 1 : p(2) + 1
                        for ik = 1 : p(1) + 1
                            N0(1, l) = N01{1}(i, ik, 1) * N01{2}(j, jk, 1) * N01{3}(k, kk, 1);
                            N1(1, l) = N01{1}(i, ik, 2) * N01{2}(j, jk, 1) * N01{3}(k, kk, 1);
                            N1(2, l) = N01{1}(i, ik, 1) * N01{2}(j, jk, 2) * N01{3}(k, kk, 1);
                            N1(3, l) = N01{1}(i, ik, 1) * N01{2}(j, jk, 1) * N01{3}(k, kk, 2);
                            l = l + 1;
                        end
                    end
                end
                [~, R1] = Rationalize(reshape(Weights, 1, []), N0, N1);
                
                ge = R1 * reshape(de, 3, [])';
                dxdxi = R1 * reshape(xyze, [], nen)';
                
                ge = dxdxi ^ (-1) * ge;
                g(:, :, i, j, k) = ge;
            end
        end
    end
end


end