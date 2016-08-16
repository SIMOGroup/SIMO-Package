function EvalPts = BsplineEval(KntVects, CtrlPts, ParaPts)
% EvalPts = BsplineEval(KntVects, CtrlPts, ParaPts)
% ------------------------------------------------------------------
% Evalutate parametric points
% -----------------------------------------------------------
% Input:
%       KntVects: knot vectors in three directions
%       CtrlPts: control_points & field
%       ParaPts: parametric points in three directions
% -----------------------------------------------------------
% Output:
%       EvalPts: matrix of interpolated points
% -----------------------------------------------------------

assert(iscell(KntVects), 'Knot vector(s) must be in cell format');
assert(iscell(ParaPts), 'Parameter points must be in cell format');

dim = numel(KntVects);
NCtrlPts = size(CtrlPts);

p = zeros(1, dim);
idx = cell(1, dim);
N0 = cell(1, dim);
NPts = cellfun(@numel, ParaPts);

for i = 1 : dim
    p(i) = numel(KntVects{i}) - NCtrlPts(i + 1) - 1;
    idx{i} = FindSpan(NCtrlPts(i + 1), p(i), ParaPts{i}, KntVects{i});
    N0{i} = BasisFuns(idx{i}, ParaPts{i}, p(i), KntVects{i});
end

if dim == 1 % curv
    EvalPts = zeros(NCtrlPts(1), NPts(1));
    for i = 1 : p(1) + 1
        EvalPts = EvalPts + bsxfun(@times, N0{1}(:, i)', CtrlPts(:, idx{1} - p + i - 1));
    end
elseif dim == 2 % surf
    EvalPts = zeros(NCtrlPts(1), NPts(2), NPts(1));
    xi_ind = idx{1} - (p(1) + 1);
    for j = 1 : p(2) + 1
        temp = zeros(NCtrlPts(1), NPts(1), NPts(2));
        eta_ind = idx{2} - (p(2) + 1) + j;
        for i = 1 : p(1) + 1
            temp = temp + bsxfun(@times, N0{1}(:, i)',...
                CtrlPts(:, xi_ind + i, eta_ind));
        end
        temp = permute(temp, [1 3 2]);
        EvalPts = EvalPts + bsxfun(@times, N0{2}(:, j)', temp);
    end
    EvalPts = permute(EvalPts, [1 3 2]);
else % volu
    EvalPts = zeros(NCtrlPts(1), NPts(3), NPts(1), NPts(2));
    x_idx = idx{1} - (p(1) + 1);
    for k = 1 : p(3) + 1
        temp2 = zeros(NCtrlPts(1), NPts(2), NPts(3), NPts(1));
        z_idx = idx{3} - (p(3) + 1) + k;
        for j = 1 : p(2) + 1
            temp1 = zeros(NCtrlPts(1), NPts(1), NPts(2), NPts(3));
            e_idx = idx{2} - (p(2) + 1) + j;
            for i = 1 : p(1) + 1
                temp1 = temp1 + bsxfun(@times, N0{1}(:, i)', CtrlPts(:, x_idx + i, e_idx, z_idx));
            end
            temp1 = permute(temp1, [1 3 4 2]);
            temp2 = temp2 + bsxfun(@times, N0{2}(:, j)', temp1);
        end
        temp2 = permute(temp2, [1 3 4 2]);
        EvalPts = EvalPts + bsxfun(@times, N0{3}(:, k)', temp2);
    end
    EvalPts = permute(EvalPts, [1 3 4 2]);
end

end