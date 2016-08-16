function ONURBS = NURBSRepara(INURBS, Dir, StartKnt, EndKnt)
% function ONURBS = NURBSRepara(INURBS, Dir, StartKnt, EndKnt)
% ------------------------------------------------------------------
% Reparameterize the knot vectors
% ------------------------------------------------------------------
% Input:
%       INURBS: input NURBS structure
%       Dir: direction to reparameter
%           (xi--Dir = 1, eta--Dir = 2, zeta--Dir = 3)
%       StartKnt: starting knot value
%       EndKnt: ending knot value
% ------------------------------------------------------------------
% Output:
%       ONURBS: output NURBS
% ------------------------------------------------------------------

p = INURBS.Order(Dir);
KntVect = INURBS.KntVect{Dir};
a = KntVect(p + 1);
b = KntVect(end - p);
if (StartKnt == a && EndKnt == b)
    ONURBS = INURBS;
else
    assert(StartKnt < EndKnt, 'StartKnot must be less than EndKnot')
    
    KntVect = KntVect - a;
    KntVect = KntVect * (EndKnt - StartKnt) / (b - a);
    KntVect = KntVect + StartKnt;
    
    finKntVect = INURBS.KntVect;
    finKntVect{Dir} = KntVect;
    
    ONURBS = CreateNURBS(finKntVect, INURBS.CtrlPts4D);
end
end