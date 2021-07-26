function [nodes, noderange] = nodes_unity_roots_d_c_r(d, c, r)
% Return d+1 nodes within the disk centered at c with radius r.

    params = [c r];
    [nodes, noderange] = nodes_unity_roots_d_cr(d, params);
end
