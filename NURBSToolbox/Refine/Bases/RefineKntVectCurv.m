function [OKntVect, OCtrlPts] = ...
    RefineKntVectCurv(p, IKntVect, ICtrlPts, InsKnts)
% function [OKntVect, OCtrlPts] = ...
%     RefineKntVectCurv(p, IKntVect, ICtrlPts, InsKnts)
% --------------------------------------------------------------
% Refine knot vector curve 
%---------------------------------------------------------------
% Input:
%       p: order of the basis functions
%       IKntVect: input knot vector
%       ICtrlPts: input control points
%       InsKnts: inserted knots
% -------------------------------------------------------------
% Output:
%       OKntVect: output knot vector
%       OCtrlPts: output control points
% -------------------------------------------------------------
% Based on Algorithm A5.4 [The NURBS BOOK, p.164]
% -------------------------------------------------------------

NCtrlPts = size(ICtrlPts, 2); % number of input control points
InsKnts  = sort(InsKnts);
NInsKnts = numel(InsKnts); % number of inserted knots
OCtrlPts = zeros(size(ICtrlPts, 1),...
    NCtrlPts + NInsKnts);
OKntVect = zeros(1, numel(IKntVect) + NInsKnts);
a = FindSpan(NCtrlPts, p, InsKnts(1), IKntVect);
b = FindSpan(NCtrlPts, p, InsKnts(NInsKnts),...
    IKntVect);
OCtrlPts(:, 1 : a - p) = ICtrlPts(:, 1 : a - p);

OCtrlPts(:, b + NInsKnts : NCtrlPts + NInsKnts) =...
    ICtrlPts(:, b : NCtrlPts);

OKntVect(1 : a) = IKntVect(1 : a);

OKntVect(b + p + 1 + NInsKnts : numel(IKntVect) + NInsKnts) =...
    IKntVect(b + p + 1 : numel(IKntVect));

i = b + p; k = b + p + NInsKnts;
for j = NInsKnts : -1 : 1
    while InsKnts(j) <= IKntVect(i) && i > a
        OCtrlPts(:, k - p - 1) = ICtrlPts(:, i - p - 1);
        OKntVect(k) = IKntVect(i); k = k - 1; i = i - 1;
    end
    OCtrlPts(:, k - p - 1) = OCtrlPts(:, k - p);
    for l = 1 : p
        ind = k - p + l;
        alfa = OKntVect(k + l) - InsKnts(j);
        if abs(alfa) == 0
            OCtrlPts(:, ind - 1) = OCtrlPts(:, ind);
        else
            alfa = alfa / (OKntVect(k + l) - ...
                IKntVect(i - p + l));
            OCtrlPts(:, ind - 1) = alfa *...
                OCtrlPts(:, ind - 1) +...
                (1 - alfa) * OCtrlPts(:, ind);
        end
    end
    OKntVect(k) = InsKnts(j); k = k - 1;
end
end