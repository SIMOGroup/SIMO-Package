function Conn = findConn(NURBS, Knt, dir)
% function Conn = findConn(NURBS, Knt, dir)
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

mcp = NURBS.NCtrlPts(1);
ncp = NURBS.NCtrlPts(2);
k = FindSpan(NURBS.NCtrlPts(1), NURBS.Order(dir), Knt, NURBS.KntVect{dir});
FDofs = sub2ind([mcp, ncp], (k - NURBS.Order(dir) - 1) * ones(1, ncp), 1 : ncp)';
SDofs = sub2ind([mcp, ncp], (k - NURBS.Order(dir)) * ones(1, ncp), 1 : ncp)';
TDofs = sub2ind([mcp, ncp], (k - NURBS.Order(dir) + 1) * ones(1, ncp), 1 : ncp)';

Conn = [FDofs SDofs; SDofs TDofs];
end