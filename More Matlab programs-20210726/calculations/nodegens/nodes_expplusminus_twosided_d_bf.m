function [nodes, noderange] = nodes_expplusminus_twosided_d_bf(d, params)
% Return d+1 nodes corresponding to powers of b from b^-d to b^0, times
% +/- f.

    b = params(1);
    f = params(2);
    nodes = [ ...
        (sym(b).^(d-1-floor(d/2):-2:-floor(d/2)) * -f) ...
        (sym(b).^(-floor(d/2):2:d-floor(d/2)) * f) ...
    ];
    noderange = nodes;
end
