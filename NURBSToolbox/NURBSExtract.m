function ONURBS = NURBSExtract(INURBS, Dir, val)
% function ONURBS = NURBSExtract(INURBS, Dir, val)
% ------------------------------------------------------------------
% Extract lower dimensional NURBS object.
%-------------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       Dir: direction along which to extract
%       (xi: dir = 1, eta: dir = 2, zeta: dir = 3)
%       val: parametric value from which to exact
%-------------------------------------------------------------------
% Output:
%       ONURBS: output NURBS structure
%-------------------------------------------------------------------

% if 'Dir'th dimension of CtrlPts4D matrix is not the last
% dimension, then permute it until it lies at last dimension
if Dir ~= INURBS.Dim
    Dirs = 1 : INURBS.Dim + 1;
    Dirs(Dir + 1) = [];
    Dirs = [Dirs, Dir + 1];
    temp = permute(INURBS.CtrlPts4D, Dirs);
else
    temp = INURBS.CtrlPts4D;
end
dim = size(temp);
temp = reshape(temp, [], dim(end));

CtrlPts = CurvPntByCornerCut(INURBS.NCtrlPts(Dir),...
    INURBS.Order(Dir), INURBS.KntVect{Dir}, temp, val);

CtrlPts = reshape(CtrlPts, dim(1 : end - 1));

KntVect = INURBS.KntVect;
KntVect{Dir} = [];
KntVect = KntVect(~cellfun(@isempty, KntVect));
ONURBS = CreateNURBS(KntVect, CtrlPts);
end