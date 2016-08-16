function [C, Stress] = StressEval(NURBS, d, ParaPts, E, nu, lab)
% [C, Stress] = StressEval(NURBS, d, ParaPts, E, nu, lab)
D = getElastMat(E, nu, lab);
F1 = reshape(d, NURBS.NNP, [])';
F2 = reshape(reshape(d, NURBS.NNP, [])', [size(F1, 1), NURBS.NCtrlPts]);
g = GradEval(NURBS.KntVect, NURBS.CtrlPts4D, ParaPts, F2);
NPts = cellfun(@numel, ParaPts);
Cw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, ParaPts);
if NURBS.Dim == 2
    Stress = zeros([3, NPts]);
    for i = 1 : NPts(1)
        for j = 1 : NPts(2)
            Strain(1) = g(1, 1, i, j); % \epsilon_{xx}
            Strain(2) = g(2, 2, i, j); % \epsilon_{yy}
            Strain(3) = g(2, 1, i, j) + g(1, 2, i, j); % 2\epsilon_{xy}
            Stress(:, i, j) = D * Strain';
        end
    end
    Weights = Cw(4, :, :);
    C = bsxfun(@rdivide, Cw(1 : 3, :, :), Weights); 
end
end