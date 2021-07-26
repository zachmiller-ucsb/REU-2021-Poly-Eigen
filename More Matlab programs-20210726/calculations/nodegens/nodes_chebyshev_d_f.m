function [nodes, noderange] = nodes_chebyshev_d_f(d, params)
% Return the d+1 Chebyshev nodes of the first kind for degree d, times a
% scaling factor f.

    f = params(1);
    nodes = zeros(1,d+1);
    for j = 1:d+1
        nodes(j) = cos((2*j-1)*pi/(2*(d+1))) * f;
    end
    noderange = [-1 1] .* f;
end
