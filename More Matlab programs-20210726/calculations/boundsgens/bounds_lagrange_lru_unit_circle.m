function [u, l] = bounds_lagrange_lru_unit_circle(d, ~, ~, ~, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Using Theorem 4.6 for bounds from document 2, only for the roots of
% unity for the unit circle.

    u = (d+1);
    l = 1/u;
end
