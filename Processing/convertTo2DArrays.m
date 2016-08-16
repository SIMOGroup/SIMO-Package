function [CtrlPts, Weights] = convertTo2DArrays(NURBS)

CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';
Weights = reshape(NURBS.Weights, 1, []);
end
