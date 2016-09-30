function [OKntVect, OCtrlPts] = DegreeElevateCurv(p, IKntVect, ICtrlPts, t)
% function [OKntVect, OCtrlPts] = DegreeElevateCurv(p, IKntVect, ICtrlPts, t)
% ---------------------------------------------------------------
% Degree elevate a curve t times
% ---------------------------------------------------------------
% Input:
%       p: order (degree) of basis functions
%       IKntVect: input knot vector
%       ICtrlPts: input control points
%       t: raise the B-spline degree t times
% ---------------------------------------------------------------
% Output:
%       OKntVect: output knot vector
%       OCtrlPts: output control points
% ---------------------------------------------------------------
% Based on Algorithm A5.9 [The NURBS BOOK, p.206].
% ---------------------------------------------------------------

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

dim = size(ICtrlPts, 1);
% coefficients for degree elevating the B\'ezier segment
BezAlfs =  zeros(p + t + 1, p + 1);
% (p + t)th-degree B\'ezier control points of the current segment
EBPts = zeros(dim, p + t + 1);
% left most control points of the next B\'ezier segment
nextBPts = zeros(dim, p - 1);
% knot insertion
alfs = zeros(p, 1);

m = numel(IKntVect); % number of knots
ph = p + t;
ph2 = floor(ph / 2);

uqKntVect = IKntVect([true; diff(IKntVect(:)) > 0]);
OKntVect = zeros(1, m + numel(uqKntVect) * t);

% Compute B\'ezier degree elevation coefficients
BezAlfs(1, 1) = 1;
BezAlfs(ph + 1, p + 1) = 1;
for i = 2 : ph2 + 1
    inv = 1 / bin(ph, i - 1);
    mpi = min(p, i - 1);
    for j = max(0, i - t - 1) + 1 : mpi + 1
        BezAlfs(i, j) = inv * bin(p, j - 1) * bin(t, i - j);
    end
end
for i = ph2 + 2 : ph
    mpi = min(p, i - 1);
    for j = max(0, i - t - 1) + 1 : mpi + 1
        BezAlfs(i, j) = BezAlfs(ph - i + 2, p - j + 2);
    end
end
mh = ph; 
kInd = ph + 2; 
r = -1; 
a = p + 1; 
b = p + 2; 
cInd = 1;
Knta = IKntVect(1);
OCtrlPts(:, 1) = ICtrlPts(:, 1);
OKntVect(1 : ph + 1) = Knta; 
% Initialize first B\'ezier seg
bPts = ICtrlPts(:, 1 : p + 1); 
while b < m % Big loop thru knot vector
    i = b;
    while b < m && IKntVect(b) == IKntVect(b + 1)
        b = b + 1;
    end
    mul = b - i + 1; 
    mh = mh + mul + t;
    Kntb = IKntVect(b); 
    oldr = r; 
    r = p - mul;
    % Insert knot knot(b) r times
    if oldr > 0
        lbz = floor((oldr + 2) / 2);
    else
        lbz = 1;
    end
    if r > 0
        rbz = ph - floor((r + 1) / 2);
    else
        rbz = ph;
    end
    if r > 0
        % Insert knot to get B\'ezier segment
        numer = Kntb - Knta;
        for k = p : -1 : mul + 1
            alfs(k - mul) = numer / (IKntVect(a + k) - Knta);
        end        
        for j = 1 : r
            s = mul + j;
            for k = p + 1 : -1 : s + 1
                bPts(:, k) = alfs(k - s) * bPts(:, k)...
                    + (1 - alfs(k - s)) * bPts(:, k - 1);
            end
            nextBPts(:, r - j + 1) = bPts(:, p + 1);
        end
    end % End of ``insert knot''    
    for i = lbz + 1 : ph + 1% Degree elevate B\'ezier
        % Only points lbz,...,ph are used below
        EBPts(:, i) = 0;
        mpi = min(p, i - 1);
        for j = max(0, i - t - 1) + 1 : mpi + 1
            EBPts(:, i) = EBPts(:, i) + BezAlfs(i, j) *...
                bPts(:, j);
        end
    end % End of degree elevating B\'ezier    
    if oldr > 1
        % Must remove knot ``knot = iniKntVect(a) oldr times''
        first = kInd - 2; 
        last = kInd; 
        den = Kntb - Knta;
        bet = floor((Kntb - OKntVect(kInd - 1)) / den);
        for tr = 1 : oldr - 1
            % Knot removal loop
            i = first; 
            j = last; 
            kj = j - kInd + 1;
            while j - i > tr
                if i - 1 < cInd
                    alf = (Kntb - OKntVect(i))/...
                        (Knta - OKntVect(i));
                    OCtrlPts(:, i) = alf *...
                        OCtrlPts(:, i) + (1 - alf) *....
                        OCtrlPts(:, i - 1);
                end
                if j - 1 >= lbz
                    if j - tr <= kInd - ph + oldr
                        gam = (Kntb - OKntVect(j - tr))...
                            / den;
                        EBPts(:, kj + 1) = gam *...
                            EBPts(:, kj + 1)+ (1 - gam) *...
                            EBPts(:, kj + 2);
                    else
                        EBPts(:, kj + 1) = bet *...
                            EBPts(:, kj + 1) + (1 - bet) *...
                            EBPts(:, kj + 2);
                    end
                end
                i = i + 1; 
                j = j - 1; 
                kj = kj - 1;
            end
            first = first - 1; 
            last = last + 1;
        end
    end % End of removing knot ``knot = iniKntVect(a)''
    if a ~= p + 1% Load the knot ``Knota''
        for i = 1 : ph - oldr
            OKntVect(kInd) = Knta; 
            kInd = kInd + 1;
        end
    end
    for j = lbz : rbz % Load control points into ``finCtrlPts''
        OCtrlPts(:, cInd + 1) = EBPts(:, j + 1);
        cInd = cInd + 1;
    end
    if b < m
        % Set up for next pass thru loop
        bPts(:, 1 : r) = nextBPts(:, 1 : r);
        bPts(:, r + 1 : p + 1) = ICtrlPts(:, b - p + r : b);
        a = b; 
        b = b + 1; 
        Knta = Kntb;
    else
        OKntVect(kInd : kInd + ph) = Kntb; 
    end
end % End while loop b < m
end

function b = bin(n, k)
% Computes the binomial coefficient.
% Author: John Burkardt
mn = min(k, n - k);
if (mn < 0) 
    b = 0;
elseif mn == 0 
    b = 1;
else
    mx = max(k, n - k); 
    b = mx + 1;
    for i = 2 : mn
        b = (b * (mx + i)) / i; 
    end
end
end
