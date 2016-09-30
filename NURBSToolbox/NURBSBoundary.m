function Bdry = NURBSBoundary(NURBS, iSide)
% function Bdry = NURBSBoundary(NURBS, iSide)
% --------------------------------------------------------------
% Extract the boundary of a NURBS structure
% --------------------------------------------------------------
% For 1D case:
% 
%         iSide = 1 o------------o iSide = 2     +----> u
%          (u = 0)                  (u = 1)
% --------------------------------------------------------------
% For 2D case:
%                 iSide = 4
%                  (v = 1)
%              o------------o                   v
%              |            |                   ^
%              |            |                   |
%   iSide = 1  |            | iSide = 2         |
%    (u = 0)   |            |  (u = 1)          +----> u
%              |            |
%              o------------o
%                 iSide = 3
%                  (v = 0)
% --------------------------------------------------------------
% For 3D case:
% 
%                    o------------o
%                   /| iSide = 6 /|
%                  / | (w = 1)  / |          w
%                 o------------o  iSide = 2  ^  v
%                 |  iSide = 4 |  |(u = 1)   | /
%       iSide = 1 |  | (v = 1) |  |          |/
%        (u = 0)  |  o-------- | -o          +----> u
%                 | /  iSide = 3 / 
%                 |/    (v = 0)|/
%                 o------------o
%                    iSide = 5   
%                     (w = 0)
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

Idx = floor((iSide + 1) / 2);
if mod(iSide, 2) == 0
    Bdry = NURBSExtract(NURBS, Idx, 1);
else
    Bdry = NURBSExtract(NURBS, Idx, 0);
end
end