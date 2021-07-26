function [u, l] = bounds_lagrange(d, nodes, ~, ~, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Using Theorem 4.1 from document 2.

    x_h = max(abs(nodes));
    x_H = max(1, x_h);
    GXd = stdgautschibound(d, nodes);
    u = (d+1)^3 * x_H^d * GXd;
    l = 1 / u;
end
