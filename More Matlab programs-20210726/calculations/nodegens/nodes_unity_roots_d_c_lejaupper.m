function [nodes, noderange] = nodes_unity_roots_d_c_lejaupper(d, params)
% Return d+1 nodes on the circle centered at c with radius based on
% Leja's upper bound.

    c = params(1);
    radius = abs(c) * (1 - (1 / 2 / (d+1)^2)^(1/d));
    [nodes, noderange] = nodes_unity_roots_d_c_r(d, c, radius);
end
