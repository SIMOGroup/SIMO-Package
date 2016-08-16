function solu = Assignment11ExactSolution(E)
% function solu = Assignement11ExactSolution(E)

solu.displ = @(x) 1 / E(x) * (56.25 - 6.25 * (1 + x) .^ 2 - 7.5 * log((1 + x) / 3));
solu.sigma = @(x) -12.5 - 12.5 * x - 2.5 ./ (1/3 + 1/3 * x);
end