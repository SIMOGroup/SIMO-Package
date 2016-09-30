function s = FindMult(Knt, KntVect)
% s = FindMult(Knt, KntVect)
% -------------------------------------------------------------------------
% Find the multiplicity of a knot value
% -------------------------------------------------------------------------
% Input:
%       Knt: knot value
%       KntVect: knot vector
% -------------------------------------------------------------------------
% Output:
%       s: multiplicity
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

s = 0;
for t = 1 : numel(KntVect)
    if Knt == KntVect(t)
        s = s + 1;
    end
end
end
