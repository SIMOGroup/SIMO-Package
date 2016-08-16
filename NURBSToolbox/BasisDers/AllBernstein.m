function B = AllBernstein(p, xi)
% function B = AllBernstein(p, xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute all pth-order Bernstein basis functions.
%--------------------------------------------------------------
% Input:
%      p: order of polynominal
%      xi: parametric points    
%--------------------------------------------------------------
% Output:
%      B: bernstein basis functions  
%--------------------------------------------------------------
% Based on Algorithm A1.3 [The NURBS BOOK, p.20]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = p + 1;
B = zeros(numel(xi), n);
Bi = zeros(n, 1);
for jj = 1 : numel(xi)
    Bi(1) = 1;
    u1 = 1 - xi(jj);
    for j = 2 : n
        saved = 0;
        for k = 1 : j - 1
            temp = Bi(k);
            Bi(k) = saved + u1 * temp;
            saved = xi(jj) * temp;
        end
        Bi(j) = saved;
    end
    B (jj, :) = Bi;
end
end