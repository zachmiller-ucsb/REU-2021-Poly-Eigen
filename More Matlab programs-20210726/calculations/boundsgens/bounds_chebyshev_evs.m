function [u, l] = bounds_chebyshev_evs(d, ~, ~, evs, L)
% Calculate a Chebyshev specific bounds for the ratio of condition numbers.

    u = ((d+1) * stdonenorm(L))/ (1 + min(abs(evs)));
    l = (1+ min(abs(evs))) / ((d+1) * stdonenorm(inv(L)));
end
