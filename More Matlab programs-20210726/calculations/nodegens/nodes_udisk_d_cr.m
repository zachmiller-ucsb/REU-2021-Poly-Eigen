function [nodes, noderange] = nodes_udisk_d_cr(d, params)
% Return d+1 nodes within the disk centered at c with radius r.

    c = params(1);
    r = params(2);
    noderange = [-1-1i 1+1i] / sqrt(2) * r + c;
    nodes = stdrandumtx(1, d+1, noderange);
end
