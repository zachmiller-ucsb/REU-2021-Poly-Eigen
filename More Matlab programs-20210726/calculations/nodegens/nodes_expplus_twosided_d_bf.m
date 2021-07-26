function [nodes, noderange] = nodes_expplus_twosided_d_bf(d, params)
% Return d+1 nodes corresponding to powers of b from b^-d to b^0, times
% +/- f.

    b = params(1);
    f = params(2);
    nodes = [ ...
        (sym(b).^(ceil(d/2)-1:-1:0) * -f) ...
        (sym(b).^(0:floor(d/2)) * f) ...
    ];
    noderange = nodes;
end
