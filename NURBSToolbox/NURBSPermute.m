function tvol = NURBSPermute(vol, ord)
% 
% NRBPERMUTE: Rearrange the directions of a NURBS volume or surface.
% 
% Calling Sequence:
% 
%   tvol = nrbpermute(vol,order)
%
% INPUT:
% 
%   vol		: NURBS volume or surface, see nrbmak.
%   order   : the order to rearrange the directions of the NURBS entity.
%
% OUTPUT:
% 
%   tvol	: NURBS volume or surface with rearranged directions.
% 
% Description:
% 
%   Utility function that rearranges the directions of a NURBS volume or
%   surface. For surfaces, nrbpermute(srf,[1 2]) is the same as
%   nrbtransp(srf). NURBS curves cannot be rearranged.
%
% Example:
%
%    nrbpermute (vol, [1 3 2])
%
%    Copyright (C) 2013 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if (~iscell(vol.KntVect))
  error('A NURBS curve cannot be rearranged.');
end

tvol = CreateNURBS(vol.KntVect(ord), permute(vol.CtrlPts4D, [1, ord + 1]));
end