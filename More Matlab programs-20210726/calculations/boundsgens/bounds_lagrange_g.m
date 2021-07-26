function [u, l] = bounds_lagrange_g(d, nodes, ~, ~, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Using Theorem 4.2 for bounds from document 2.
% Using Lemma 21 for L_L and L_H expressions.

    x_h = max(abs(nodes));
    x_H = max(1, x_h);

    GXd = stdgautschibound(d, nodes);

    L_H = (d+1) * GXd;
    L_L = 1 / ((d+1) * x_H^d);

    u = (d+1) * x_H^d * L_H;
    l = L_L / (d+1) / GXd;
end
