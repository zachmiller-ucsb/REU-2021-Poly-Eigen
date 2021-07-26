function [nodes, noderange] = nodes_expminus_twosided_d_bf(d, params)
% Return d+1 nodes corresponding to powers of b from b^-d to b^0, times
% +/- f.

    b = params(1);
    f = params(2);
    nodes = [ ...
        (sym(b).^(0:-1:-ceil(d/2)+1) * -f) ...
        (sym(b).^(-floor(d/2):0) * f) ...
    ];
    noderange = nodes;
end
