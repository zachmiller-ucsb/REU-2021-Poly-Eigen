function [nodes, noderange] = nodes_expplus_d_bf(d, params)
% Return d+1 nodes corresponding to powers of b from b^0 to b^d, times f.

    b = params(1);
    f = params(2);
    nodes = sym(b).^(0:d) * f;
    noderange = nodes;
end
