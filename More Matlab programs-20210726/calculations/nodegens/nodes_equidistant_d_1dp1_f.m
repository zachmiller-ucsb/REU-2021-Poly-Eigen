function [nodes, noderange] = nodes_d_1dp1_f(d, params)
% Return d+1 equidistant nodes between 1 and d+1, times f.

    f = params(1);
    [nodes, noderange] = nodes_equidistant_d_abf(d, [1 d+1 f]);
end
