function [u, l] = bounds_chebyshev_cond_2(d, nodes, ~, ~, ~)
% Calculate a Chebyshev specific bounds for the ratio of condition numbers.

    V = stdvander(sym(nodes));
    u = (d+1)^3 * sqrt(2/pi) * stdcond(V);
    l = 1/u;
end
