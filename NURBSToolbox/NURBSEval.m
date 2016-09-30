function varargout = NURBSEval(NURBS, ParaPts, varargin)
% function varargout = NURBSEval(NURBS, ParaPts, varargin)
% ------------------------------------------------------------------
% Interpolate parameter points and field
% ------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       ParaPts: parameter points
%       varargin: field variable (optional)
% ------------------------------------------------------------------
% Output:
%       varargout: interpolated points and field (optional)
% ------------------------------------------------------------------

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

assert(iscell(ParaPts), 'ParaPts must be stored in cell format');
if (nargin == 2)
    EvalCw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, ParaPts);
    if NURBS.Dim == 1
        Weights = EvalCw(4, :);
        C = bsxfun(@rdivide, EvalCw(1 : 3, :), Weights);
    elseif NURBS.Dim == 2
        Weights = EvalCw(4, :, :);
        C = bsxfun(@rdivide, EvalCw(1 : 3, :, :), Weights);
    elseif NURBS.Dim == 3
        Weights = EvalCw(4, :, :, :);
        C = bsxfun(@rdivide, EvalCw(1 : 3, :, :, :), Weights);
    end
else
    Field = reshape(varargin{1}, NURBS.NNP, [])';
    Field = reshape(Field, [size(Field, 1), NURBS.NCtrlPts]);
    Fw = bsxfun(@times, Field, NURBS.Weights);
    CFw = cat(1, NURBS.CtrlPts4D, Fw);
    EvalCFw = BsplineEval(NURBS.KntVect, CFw, ParaPts);
    if NURBS.Dim == 1
        Weights = EvalCFw(4, :);
        C = bsxfun(@rdivide, EvalCFw(1 : 3, :), Weights);
        F = bsxfun(@rdivide, EvalCFw(5 : end, :), Weights);
    elseif NURBS.Dim == 2
        Weights = EvalCFw(4, :, :);
        C = bsxfun(@rdivide, EvalCFw(1 : 3, :, :), Weights);
        F = bsxfun(@rdivide, EvalCFw(5 : end, :, :), Weights);
    elseif NURBS.Dim == 3
        Weights = EvalCFw(4, :, :, :);
        C = bsxfun(@rdivide, EvalCFw(1 : 3, :, :, :), Weights);
        F = bsxfun(@rdivide, EvalCFw(5 : end, :, :, :), Weights);
    end
end

if (nargout == 1)
    varargout{1} = C;
elseif (nargout == 2)
    if (nargin == 3)
        varargout{1} = C;
        varargout{2} = squeeze(F);
    else
        error('You must input vector of temperature or displacement field');
    end
end
end