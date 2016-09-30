function varargout = convertToTripletStorage(Mesh, KVals, varargin)
% varargout = convertToTripletStorage(Mesh, KVals, varargin)
% -----------------------------------------------------------------------
% Input:
%       Mesh: mesh structure
%       KVals: matrix stores local stiffness matrices
%           (size(KVals) = [(NEN * DOF) ^ 2, NEL])
%       varargin: matrix stores local body forces (optional)
%           (size(FVals) = [NEN * DOF, NEL])
% ------------------------------------------------------------------------
% Output:
%       varargout: triplet sparse storage (Rows, Cols, Vals) 
%           varargout{1} = Rows;
%           varargout{2} = Cols;
%           varargout{3} = Vals;
%           varargout{4} = F;
%           ( K = sparse(Rows, Cols, Vals) )
% ------------------------------------------------------------------------

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

J = repmat(1 : Mesh.NEN * Mesh.Dof, Mesh.NEN * Mesh.Dof, 1);
I = J';
if Mesh.Dof == 1
    ElConn = Mesh.El;
elseif Mesh.Dof == 2
    ElConn = [Mesh.El, Mesh.El + Mesh.NDof / Mesh.Dof];
elseif Mesh.Dof == 3
    ElConn = [Mesh.El, Mesh.El + Mesh.NDof / Mesh.Dof, Mesh.El + Mesh.NDof / Mesh.Dof * 2];
end
ii = ElConn(:, I(:))';
jj = ElConn(:, J(:))';

varargout{1} = ii(:); % Rows
varargout{2} = jj(:); % Cols
varargout{3} = KVals(:); % Vals
if (nargin == 3)
    kk = ElConn';
    FVals = varargin{1};
    varargout{4} = accumarray(kk(:), FVals(:)); % F
end
end