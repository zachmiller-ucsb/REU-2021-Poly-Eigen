function [nodes, noderange] = nodes_equidistant_d_pmd_f(d, params)
% Return d+1 equidistant nodes between -d and d, times f.

    f = params(1);
    [nodes, noderange] = nodes_equidistant_d_abf(d, [-d d f]);
end
