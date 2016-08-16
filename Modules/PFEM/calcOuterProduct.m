function T = calcOuterProduct(varargin)

T = varargin{1};
for i = 2 : numel(varargin)
    sa = size(T);
    sb = size(varargin{i});
    S = permute(varargin{i}, circshift(1 : (ndims(T) + ndims(T)), [0, ndims(T)]));
    tmp = bsxfun(@times, T, S);
    N = ndims(tmp);
    switch N
        case 4 % one element
            if i == numel(varargin)
                tmp = permute(tmp, [2, 4, 1, 3]);
            else
                tmp = permute(tmp, [1, 3, 2, 4]);
            end
        case 6 % more than one element
            if i == numel(varargin)
                tmp = permute(tmp, [2, 5, 1, 4, 3, 6]);
            else
                tmp = permute(tmp, [1, 4, 2, 5, 3, 6]);
            end
    end
    if i == numel(varargin)
        T = reshape(tmp, 1, sa(2) * sb(2), sa(1) * sb(1), []);
    else
        T = reshape(tmp, sa(1) * sb(1), sa(2) * sb(2), []);
    end
end
end