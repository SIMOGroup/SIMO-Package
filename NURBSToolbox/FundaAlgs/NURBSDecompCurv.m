function [nb, Qw] = NURBSDecompCurv(n, p, KntVect, Pw)
% [nb, Qw] = NURBSDecompCurv(n, p, KntVect, Pw)
% -------------------------------------------------------------------------
% Decompose curve into Bezier segments
% ---------------------------------------------------------------
% Input: 
%       n: number of control points.
%       p: order (degree) of basis functions
%       KntVect: knot vector
%       Pw: control points
% ---------------------------------------------------------------
% Ouput: 
%       nb: number of segments.
%       Qw: control points of bezier segments.
% ---------------------------------------------------------------
% Based on Algorithm A5.6 from ``The NURBS BOOK'' pg173.
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

m = n + p;
a = p;
b = p + 1;
nb = 0;

Qw{nb + 1} = Pw(:, 1 : p + 1);

while b < m
    i = b;
    while b < m && KntVect(b + 1 + 1) == KntVect(b + 1)
        b = b + 1;
    end
    mult = b - i + 1;
    if mult < p
        numer = KntVect(b + 1) - KntVect(a + 1); % Numerator of alpha
        % compute and store alphas
        for j = p : -1 : mult + 1
            alphas(j - mult) = numer / (KntVect(a + j + 1) - KntVect(a + 1));
        end
        r = p - mult;   % insert knot r times
        for j = 1 : r
            save = r - j;
            s = mult + j;
            for k = p : -1 : s
                alpha = alphas(k - s + 1);
                Qw{nb + 1}(:, k + 1) = alpha * Qw{nb + 1}(:, k + 1) + (1 - alpha) * Qw{nb + 1}(:, k);
            end
            if b < m
                Qw{nb + 2}(:, save + 1) = Qw{nb + 1}(:, p + 1); % next segment
            end
        end
    end
    nb = nb + 1;
    if b < m
        for i = p - mult : p
            Qw{nb + 1}(:, i + 1) = Pw(:, b - p + i + 1);
        end
        a = b;
        b = b + 1;
    end
end % while
end
