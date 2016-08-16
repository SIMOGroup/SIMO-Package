function NURBS = CreateNURBS(KntVect, CtrlPts)
% NURBS = CreateNURBS(KntVects, CtrlPts)
% ------------------------------------
assert(iscell(KntVect), 'Input knots must be in cell format');
assert(numel(KntVect) >= 1);
assert(numel(KntVect) <= 3);
for i = 1 : numel(KntVect)
    assert(size(KntVect{i}, 1) == 1)
    assert(numel(KntVect{i}) >= 4, ...
        'Number of knot values must be equal or greater than 4')
end
% ------------------------------------
Dim = numel(KntVect);

if numel(KntVect) == 1
    W = CtrlPts(4, :);
    CtrlPts3D = bsxfun(@rdivide, CtrlPts(1 : 3, :), W);
elseif numel(KntVect) == 2
    W = CtrlPts(4, :, :);
    CtrlPts3D = bsxfun(@rdivide, CtrlPts(1 : 3, :, :), W);
else
    W = CtrlPts(4, :, :, :);
    CtrlPts3D = bsxfun(@rdivide, CtrlPts(1 : 3, :, :, :), W);
end

NPts = size(CtrlPts);
NCtrlPts = zeros(1, Dim);
p = zeros(1, Dim);
uqKntVect = cell(1, Dim);
for i = 1 : Dim
    NCtrlPts(i) = NPts(i + 1);
    p(i) = numel(KntVect{i}) - NCtrlPts(i) - 1;
    uqKntVect{i}  = KntVect{i}([true; diff(KntVect{i}(:)) > 0]);
end
NURBS.KntVect = KntVect;
NURBS.uqKntVect = uqKntVect; % unique knot values
NURBS.CtrlPts4D = CtrlPts; % control points in 4D space
NURBS.CtrlPts3D = CtrlPts3D; % control points projected into 3D space
NURBS.Weights = W;
NURBS.Dim = Dim;
% number of control points in each direction
NURBS.NCtrlPts = NCtrlPts;
NURBS.Order = p;
% number of total control points
NURBS.NNP = prod(NCtrlPts); %"NP" is an abbreviation for nodal points
end