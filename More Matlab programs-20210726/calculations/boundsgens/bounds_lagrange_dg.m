function [u, l] = bounds_lagrange_dg(d, nodes, ~, evs, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Using Theorem 4.2 for bounds from document 2.
% Using Lemma 21 for L_L and L_H expressions.

    [ud, ld] = bounds_lagrange_d(d, nodes, [], evs, []);
    [ug, lg] = bounds_lagrange_g(d, nodes, [], evs, []);
    u = min(ud, ug);
    l = max(ld, lg);
end
