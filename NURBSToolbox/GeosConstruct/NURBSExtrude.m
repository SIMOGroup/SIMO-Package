function ONURBS = NURBSExtrude(INURBS, DisplVect)
% ONURBS = NURBSExtrude(INURBS, DisplVect)
% --------------------------------------------------------------------
% Construct a NURBS surface (volume) by extruding a NURBS curve (surface)
% ----------------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       DisplVect: displacment vector
%       (e.g. DisplVect = [0 0 11] means the height of the extrusion
%              is 11 units along the z axis)
% ----------------------------------------------------------------------
% Ouput:
%       ONURBS: input NURBS structure
% ----------------------------------------------------------------------

% translation matrix
A = eye(4);
A(1 : end - 1, 4) = DisplVect';

a = INURBS.CtrlPts4D;
b = reshape(A * reshape(a, 4, []), size(a));

CtrlPts = cat(INURBS.Dim + 2, a, b);

ONURBS = CreateNURBS([INURBS.KntVect, [0, 0, 1, 1]], CtrlPts);
end