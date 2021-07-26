function [nodes, noderange] = nodes_equidistant_d_abf(d, params)
% Make d+1 equidistant nodes between a and b, times f.

    a = params(1);
    b = params(2);
    f = params(3);
    nodes = ((sym(b) - sym(a)) / d * (0:d) + sym(a)) * f;
    noderange = [a b] .* f;
end
