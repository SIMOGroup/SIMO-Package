% initialize
clear all; close all; clc;

% choose parameters of quadrature rule
p = 10;              % choose polynomial degree
r = 3;               % regularity in between the elements

% test rule for elements
passed = true;
for e=1:100

    % compute the quadrature rule
    [nodes, weights] = quadrule(p,r,e,0.0,e);

    % compute analytical integrals
    knots  = [0.0*ones(1,p) linspace(0.0,e,e+1) e*ones(1,p)];
    dim = length(knots)-p-1;
    target = zeros(dim,1);
    for k=1:dim
        target(k) = (knots(k+p+1)-knots(k)) / (p+1);
    end
    
    % compute basisfunctions
    B = spcol(knots,p+1,nodes);
    
    % normalized error
    err = norm(target' - weights' * B)/dim;

    % check if quadrature rule passes
    if err<10e-15
        passed = passed*true;
    else
        passed = false;
    end
end

% output result
sprintf('rule (%d,%d) passed: %d',p,r,passed)