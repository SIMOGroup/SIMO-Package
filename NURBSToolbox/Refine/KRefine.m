function ONURBS = KRefine(INURBS, NEl, Order, Continuity)
% function ONURBS = KRefine(INURBS, NEl, Order, Continuity)
% -----------------------------------------------------------------
% Refine a NURBS object (curve, surface, volume) by knot refinement
% and order (degree) elevation
% -----------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       NEl: number of elements per direction
%       Order: order of basis functions per direction
%       Continuity: continuity of basis functions per direction
% -------------------------------------------------------------
% Output:
%       ONURBS: output NURBS structure
% -------------------------------------------------------------

ONURBS = INURBS;
% order elevation
assert(all(Order - ONURBS.Order >= 0));
if all(Order)
    for Dir = 1 : ONURBS.Dim
        p = Order(Dir);
        t = p - ONURBS.Order(Dir);
        if t > 0
            ONURBS = PRefine(ONURBS, Dir, t);
        end
    end
end
% knot refinement
assert(all(Continuity >= 0));
assert(all(Continuity < Order));
assert(all(NEl > 0));
for Dir = 1 : ONURBS.Dim
    [Knts, Mlts] = knt2brk(ONURBS.KntVect{Dir});
    N = NEl(Dir);
    % compute the number of knots to be inserted
    deltaXi = diff(Knts) / N;
    Knts = Knts(1 : end - 1);
    
    Xi = repmat(Knts', 1, N);
    step = 1 : N - 1;
    Xi(:, 2 : end) = Xi(:, 2 : end) + deltaXi' * step;
    
    Mlts = Mlts(1 : end - 1);
    repsMat = zeros(numel(Mlts), N);
    repsMat(:, 1) = ONURBS.Order(Dir) - Continuity(Dir) - Mlts;
    repsMat(:, 2 : end) = ONURBS.Order(Dir) - Continuity(Dir);
    
    Knts = reshape(Xi', 1, []);
    repsVect = reshape(repsMat', 1, []);
    
    Knts = Knts(2 : end);
    repsVect = repsVect(2 : end);
    
    insKnts = Knts(repsVect > 0); % inserted knots
    repsVect = repsVect(repsVect > 0);
    
    Idx = zeros(1, sum(repsVect));
    Idx(cumsum([1 repsVect(1 : end - 1)])) = 1;
    Idx = cumsum(Idx);
    
    if ~isempty(insKnts)
        KntsMult = insKnts(Idx);
        % insert multiple knots
        ONURBS = HRefine(ONURBS, Dir, KntsMult);
    end
end
end