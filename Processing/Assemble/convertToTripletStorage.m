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