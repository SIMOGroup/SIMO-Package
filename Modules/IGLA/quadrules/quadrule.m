function [nodes, weights] = quadrule(p,r,e,a,b)
% this function computes the quadrature rule for the target splinspace 
% S^p_1 on a uniform partition with open knotvector.
% input:
%      p - polynomial degree
%      r - internal regularity in between elements
%      m - number of elements
%      a - left boundary of parametric domain (first knot)
%      b - right boundary of parametric domain (last knot)

% initialize quadrature rules
n = (p+1)*2 + (e-1)*(p-r) - p-1;   % dimension of the spline space
m = int64(ceil(n/2));              % dimension of the quadrature rule
nodes   = zeros(m,1);                % allocate space for nodes
weights = zeros(m,1);                % allocate space for weights


 % load data & initialize
if p==4 && r==0
    if mod(e,2)==0
        general_rule = sprintf('%s%d%s%d%s%d%s%d%s','quadrule_p',p,'_r',r,'/quadrule_even_r',r,'_p',p,'.mat');
    else
        general_rule = sprintf('%s%d%s%d%s%d%s%d%s','quadrule_p',p,'_r',r,'/quadrule_odd_r',r,'_p',p,'.mat');
    end
else
    if mod(m,2)==0
        general_rule = sprintf('%s%d%s%d%s%d%s%d%s','quadrule_p',p,'_r',r,'/quadrule_even_r',r,'_p',p,'.mat');
    else
        general_rule = sprintf('%s%d%s%d%s%d%s%d%s','quadrule_p',p,'_r',r,'/quadrule_odd_r',r,'_p',p,'.mat');
    end
end
load(general_rule); m1 = length(bnodes); m2 = length(inodes); m3 = length(mnodes);

% rule for even number of nodes
if 2*m==n

    % general rule with repetetive structure
    if m>2*m1

        % construct quadrature rule in the uniform parameter domain with unit spacing
        z = ceil(0.5*m);
        nodes(1:m1)   = bnodes;
        weights(1:m1) = bweights;
        left  = ceil(bnodes(m1));
        right = ceil(inodes(m2));
        ii = m1;
        while ii<z
            nodes(ii+1:ii+m2)   = inodes + left;
            weights(ii+1:ii+m2) = iweights;
            left                = left + right;
            ii                  = ii+m2;
        end
        nodes(m-z+1:m) = e - flipud(nodes(1:z));
        weights(m-z+1:m) = flipud(weights(1:z));

    % load precomputed specific rule    
    else
        specific_rule = sprintf('%s%d%s%d%s%d%s%d%s%d%s','quadrule_p',p,'_r',r,'/quadrule_e',e,'_r',r,'_p',p,'.mat');
        load(specific_rule);
    end

% rule for odd number of nodes
else

    % rule with repetetive structure
    if m>2*(m1+m3)-1

        % construct quadrature rule in the uniform parameter domain with unit spacing
        z = ceil(0.5*m);
        nodes(1:m1)   = bnodes;
        weights(1:m1) = bweights;
        left  = ceil(bnodes(end));
        right = ceil(inodes(end));
        ii = m1;
        while ii<z
            nodes(ii+1:ii+m2)   = inodes + left;
            weights(ii+1:ii+m2) = iweights;
            left                = left + right;
            ii                  = ii+m2;
        end
        nodes(z-m3+1:z)   = 0.5*e + mnodes;
        weights(z-m3+1:z) = mweights; 

        nodes(m-z+1:m) = e - flipud(nodes(1:z));
        weights(m-z+1:m) = flipud(weights(1:z));

    % load precomputed rule
    else
        specific_rule = sprintf('%s%d%s%d%s%d%s%d%s%d%s','quadrule_p',p,'_r',r,'/quadrule_e',e,'_r',r,'_p',p,'.mat');
        load(specific_rule);
    end
end
    
 
% rescale rule for elementsize h
h = (b-a) / e;
nodes = a + h.*nodes;
weights = h * weights;

end




