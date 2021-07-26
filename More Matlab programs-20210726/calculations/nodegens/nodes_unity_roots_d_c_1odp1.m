function [nodes, noderange] = nodes_unity_roots_d_c_1odp1(d, params)
% Return d+1 nodes within the disk centered at c with radius 1/(d+1).

    newparameters = [params(1) (1 / (d+1))];
    [nodes, noderange] =  nodes_unity_roots_d_cr(d, newparameters);
end
