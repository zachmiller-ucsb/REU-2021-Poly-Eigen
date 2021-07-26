function [u, l] = bounds_chebyshev_dp1_evs(d, ~, ~, evs, L)
% Calculate a Chebyshev specific bounds for the ratio of condition numbers.

    u = (d+1)^2 * max(1 ./ (1 + abs(evs)));
    l = 1 / ((d+1) * stdonenorm(inv(L)));
end
