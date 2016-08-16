function solu = CantiliverBeamExactSolu(L, D, E, I, nu, P)

solu.ux = @(x, y) P * y / (6 * E * I) .* ((6 * L - 3 * x) .* x ... 
          + (2 + nu) * (y .^ 2 - D ^ 2 / 4));
solu.uy = @(x, y) - P / (6 * E * I) * (3 * nu * y .^ 2 .* (L - x)... 
          + (4 + 5 * nu) * D ^ 2 .* x / 4 + (3 * L - x) .* x .^ 2);
end