function [u, l] = bounds_chebyshev_upperbounds(d, ~, ~, evs, L)
% Chebyshev specific upperbounds for the ratio of condition numbers.
% Note this compares the evs upper bounds to (d+1)^2.

    [u, ~] = bounds_chebyshev_evs(d, [], [], evs, L);
    l = (d+1)^2;
end
