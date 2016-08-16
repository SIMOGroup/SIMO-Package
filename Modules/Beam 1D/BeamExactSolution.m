function solu = BeamExactSolution(L, k, q, F)
% function solu = BeamExactSolution(E, I, L, k, q, F)
% Calculate exact solution for bar problem

solu.displ = @(x)F/(6*k)*(3*L*x.^2-x.^3);
%solu.moment = @(x)q/(24*k)*(x.^4 - 4*L*x.^3 + 6*L^2*x.^2); for
%distribution
end