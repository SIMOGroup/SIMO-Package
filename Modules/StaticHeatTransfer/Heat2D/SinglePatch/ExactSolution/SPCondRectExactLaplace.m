function T = SPCondRectExactLaplace(fx, x, y, a, b, n)
% T = SPCondRectExactLaplace(fx, x, y, a, b, n)
% Evaluate Laplace solution
% ---------------------------------------------------------------
% Input:
%       fx: function on boundary y = b
%       x: x value, 0 <= x <= a
%       y: y value, 0 <= y <= b
%       a: first dimension of the rectangular domain
%       b: second dimension of the rectangular domain
%       n: number of summation
% -------------------------------------------------------------
% Output:
%       T: solution of laplace equation at given values
% -------------------------------------------------------------
nArray = 1 : n;
fun = @(x) fx(x) * sin((nArray * pi * x) / a);
bn = (2 / a * integral(fun, 0, a, 'ArrayValued',true)) ./ (sinh((nArray * pi * b) / a));
tmp = zeros([size(x), n]);
for i = 1 : n
    tmp(:, :, i) = bn(i) * sin(i * pi * x/a) .* sinh(i * pi * y/a);
end
T = sum(tmp, 3);
end