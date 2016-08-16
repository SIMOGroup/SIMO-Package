function Cw = NURBSLine(P1, P2)
% Cw = NURBSLine(P1, P2)
% -----------------------------------------------------------------
% Contruct a straight line
% -----------------------------------------------------------------
% Input:
%       P1: the 1st point of the line
%       P2: the 2nd point of the line
% -----------------------------------------------------------------
% Output:
%       Cw: NURBS structure of the line
% -----------------------------------------------------------------

CtrlPts = zeros(4, 2);
CtrlPts(1 : numel(P1), 1) = P1;
CtrlPts(1 : numel(P2), 2) = P2;
CtrlPts(4, :) = 1;
KntVect = [0 0 1 1];
Cw = CreateNURBS({KntVect}, CtrlPts);
end