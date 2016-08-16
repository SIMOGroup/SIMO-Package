function s = FindMult(Knt, KntVect)
% s = FindMult(Knt, KntVect)
% -------------------------------------------------------------
% Find the multiplicity of a knot value
% -------------------------------------------------------------
% Input:
%       Knt: knot value
%       KntVect: knot vector
% -------------------------------------------------------------
% Output:
%       s: multiplicity
% -------------------------------------------------------------

s = 0;
for t = 1 : numel(KntVect)
    if Knt == KntVect(t)
        s = s + 1;
    end
end
end
