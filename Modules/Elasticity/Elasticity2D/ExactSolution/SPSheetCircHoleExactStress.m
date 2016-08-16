function S = SPSheetCircHoleExactStress(a, Tx)
% S = SPSheetCircHoleExactStress(a, Tx)
% Evaluate analytical solution of stresses for "sheet with circular hole"
% problem
% ------------------------------------------------------------------
% Input:
%       a: radius of the hole
%       Tx: uniaxial tension
% ------------------------------------------------------------------
% Output:
%       S: a structure contains 3 stress functions,
%           S.xx(x, y)
%           S.yy(x, y)
%           S.xy(x, y)
% ------------------------------------------------------------------
r = @(x, y) sqrt(x .^ 2 + y .^ 2);
theta = @(x, y) atan(y ./ x);
C2t = @(x, y) cos(2 * theta(x, y));
C4t = @(x, y) cos(4 * theta(x, y));
S2t = @(x, y) sin(2 * theta(x, y));
S4t = @(x, y) sin(4 * theta(x, y));
term1 = @(x, y) (a ./ r(x, y)) .^ 2;
term2 = @(x, y) (a ./ r(x, y)) .^ 4;
S.xx = @(x, y) Tx * (1 - term1(x, y) .* (3 / 2 * C2t(x, y) +...
    C4t(x, y)) + 3 / 2 * term2(x, y) .* C4t(x, y));
S.yy = @(x, y) Tx * (-term1(x, y) .* (1 / 2 * C2t(x, y) -...
    C4t(x, y)) - 3 / 2 * term2(x, y) .* C4t(x, y));
S.xy = @(x, y) Tx * (-term1(x, y) .* (1 / 2 * S2t(x, y) +...
    S4t(x, y)) + 3 / 2 * term2(x, y) .* S4t(x, y));
end