function [nodes, noderange] = nodes_unity_roots_d_c_sinpiodp1(d, params)
% Return d+1 nodes on the circle centered at 0 with radius sin(pi/(d+1)).

    newparams = [params(1) sin(pi/(d+1))];
    [nodes, noderange] = nodes_unity_roots_d_c_r(d, 0, newparams);
end
