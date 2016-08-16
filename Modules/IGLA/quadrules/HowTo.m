% initialize
clear; close all; clc;

% choose parameters
p = 3;              % choose polynomial degree
r = 0;              % regularity in between the elements
e = 12;             % number of elements (an even number are required for the even degree cases)
a = 0.0;            % left boundary
b = e;              % right boundary

% compute the quadrature rule
[nodes, weights] = quadrule(p,r,e,a,b);

% plot the quadrature nodes and weights
stem(nodes,weights,'o')
xlabel('x')
ylabel('w')
title('quadrature nodes versus weights')

% test rule
knots  = [a*ones(1,p) linspace(a,b,e+1) b*ones(1,p)];
dim = length(knots)-p-1;
target = zeros(dim,1);
for k=1:dim
    target(k) = (knots(k+p+1)-knots(k)) / (p+1);
end
B = spcol(knots,p+1,nodes);
err = norm(target' - weights' * B)/dim;

% output result
sprintf('rule (%d,%d) passed: %d',p,r,err<10e-15)
