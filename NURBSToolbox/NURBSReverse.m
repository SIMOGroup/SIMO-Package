function nrb = NURBSReverse(nrb, idir)
%
% NURBSReverse: Reverse the evaluation directions of a NURBS geometry.
% 
% Calling Sequence:
% 
%   rnrb = nrbreverse(nrb);
%   rnrb = nrbreverse(nrb, idir);
% 
% INPUT:
% 
%   nrb		: NURBS data structure, see nrbmak.
%   idir        : vector of directions to reverse.
%
% OUTPUT:
% 
%   rnrb	: Reversed NURBS.
% 
% Description:
% 
%   Utility function to reverse the evaluation direction of a NURBS
%   curve or surface.
%
%    Copyright (C) 2000 Mark Spink
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

if (nargin > 2)
  error('Incorrect number of input arguments');
end

if (iscell(nrb.KntVect))
  % reverse a NURBS surface or volume
  ndim = numel (nrb.KntVect);
  if (nargin == 1 || isempty (idir))
    idir = 1:ndim;
  end
  for ii = idir
    nrb.KntVect{ii} = sort(nrb.KntVect{ii}(end) - nrb.KntVect{ii});
    nrb.CtrlPts4D = flip(nrb.CtrlPts4D, ii + 1);
    nrb.CtrlPts3D = flip(nrb.CtrlPts3D, ii + 1);
  end

else
  % reverse a NURBS curve
  nrb.KntVect = sort(nrb.KntVect(end) - nrb.KntVect);
  nrb.CtrlPts4D = fliplr(nrb.CtrlPts4D);
  nrb.CtrlPts3D = fliplr(nrb.CtrlPts3D);
end

end
