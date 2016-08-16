function SPToVTK(NURBS, d, ParaPts, filename, fieldname)
% SPToVTK(NURBS, d, ParaPts, filename, fieldname)
% --------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       d: temperature or displacement field
%       ParaPts: parameter points per direction
%       filename: name of file
%       fieldname: name of field
% ------------------------------------------------------------------
% Output:
%       filename.vts
% ------------------------------------------------------------------

F1 = reshape(d, NURBS.NNP, [])';
F2 = reshape(F1, [size(F1, 1), NURBS.NCtrlPts]);
F2w = bsxfun(@times, F2, NURBS.Weights);
CFw = cat(1, NURBS.CtrlPts4D, F2w);
tmp = BsplineEval(NURBS.KntVect, CFw, ParaPts);
Weights = tmp(4, :, :, :);
C = bsxfun(@rdivide, tmp(1 : 3, :, :, :), Weights);
F = bsxfun(@rdivide, tmp(5 : end, :, :, :), Weights);
exportToVTK(C, squeeze(F), filename, fieldname);
end