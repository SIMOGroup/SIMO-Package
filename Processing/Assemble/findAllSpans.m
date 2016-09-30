function idcs = findAllSpans(n, p, KntVect, NEl)
% function idcs = findAllSpans(n, p, KntVect, NEl)
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

idcs = zeros(1, NEl); % span index
el = 1;
for i = p + 1 : n
    if (abs(KntVect(i) - KntVect(i + 1)) > sqrt(eps))
        idcs(el) = i;
        el = el + 1;
    end
end
end