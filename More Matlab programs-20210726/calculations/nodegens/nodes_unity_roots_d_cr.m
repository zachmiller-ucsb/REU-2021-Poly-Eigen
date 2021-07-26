function [nodes, noderange] = nodes_unity_roots_d_cr(d, params)
% Return d+1 nodes within the disk centered at c with radius r, given as
% paramaters.

    c = params(1);
    r = params(2);
    noderange = [-1-1i 1+1i] / sqrt(2) * r + c;
    nodes = exp((0:d) * 2i * pi / (d+1)) * r + c;
end
