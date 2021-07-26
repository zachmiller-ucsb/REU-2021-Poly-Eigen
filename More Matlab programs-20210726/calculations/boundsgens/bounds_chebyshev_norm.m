function [u, l] = bounds_chebyshev_norm(d, ~, ~, ~, L)
% Calculate a Chebyshev specific bounds for the ratio of condition numbers.

    u = (d+1) * stdonenorm(L);
    l = 1 / ((d+1) * stdonenorm(inv(L)));
end
