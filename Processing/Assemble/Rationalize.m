function varargout = Rationalize(W, N0, N1, varargin)
% function varargout = Rationalize(W, N0, N1, varargin)
% ----------------------------------------------------------------------
% Rationalize Bspline basis functions and their first derivatives
% ----------------------------------------------------------------------
% Input:
%       W: weights of control points, size(W) = [1, NEN]
%       N0: bspline basis functions, size(N0) = [1, NEN]
%       N1: first derivatives of bsplines, size(N1) = [Dim, NEN],
%           Dim is dimension of problem (Dim \neq NSD)
%       varargin (N2-optional): second derivatives of bsplines, 
%           for 1D, size(N2) = [1, NEN],
%           for 2D, size(N2) = [3, NEN],
%-----------------------------------------------------------------------
% Output:
%       varargout{1} = R0: NURBS basis functions
%       varargout{2} = R1: corresponding first derivatives
%       varargout{3} = R2: corresponding second derivatives
%           varargout{3}(1, :) = d2Rdxi2
%           varargout{3}(2, :) = d2Rdeta2
%           varargout{3}(3, :) = d2Rdxideta
% ----------------------------------------------------------------------

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

% Convert Bspline to NURBS basis functions for 1D, 2D and 3D
N0W = N0 .* W;
W0 = sum(N0W);
R0 = N0W / W0;
varargout{1} = R0;
% -------------------------------------------------------------------------
% First derivatives of NURBS basis functions for 1D, 2D and 3D
N1W = bsxfun(@times, N1, W);
W1 = sum(N1W, 2);
N1WR0W1 = N1W - bsxfun(@times, R0, W1);
R1 = N1WR0W1 / W0;
varargout{2} = R1;
% -------------------------------------------------------------------------
% Second derivatives of NURBS basis functions for 1D and 2D
if (nargin == 4)
    N2 = varargin{1};
    if(size(N2, 1) == 1)% 1D problem
        N2W = N2 .* W;
        W2 = sum(N2W);
        varargout{3} = (N2W - 2 * R1 * W1 - R0 * W2) / W0;
    elseif (size(N2, 1) == 3)% 2D problem
        N2W = bsxfun(@times, N2, W);
        W2 = sum(N2W, 2);
        varargout{3}(1 : 2, :) = (N2W(1 : 2, :) - 2 * bsxfun(@times, R1, W1) - bsxfun(@times, R0, W2(1 : 2, :))) / W0;
        varargout{3}(3, :) = (N2W(3, :) - R1(1, :) * W1(2) - R1(2, :) * W1(1) - R0 * W2(3, :)) / W0;
    end
end
end