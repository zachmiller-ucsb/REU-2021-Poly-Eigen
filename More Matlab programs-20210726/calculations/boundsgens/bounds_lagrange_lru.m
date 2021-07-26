function [u, l] = bounds_lagrange_lru(d, nodes, noderange, evs, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Using Theorem 4.2 for bounds from document 2 (improvement of previous
% case for bounds_lagrange_d for lower bound and upper bound, only
% applicable to nodes that are roots of unity).

    x_h = max(abs(nodes));
    x_H = max(1, x_h);

    % Determine the shape of the noderange.
    [~, ~, ~, radius] = stdexplainrange(noderange);
    radius = abs(radius);

    L_H_evs = max(1 ./ (1 + abs(evs).^d));
    L_L_evs = min(1 ./ max(1, abs(evs)).^d);

    u = (d+1)^2 / sin(pi / (d+1))^d * x_H^d * L_H_evs;
    l = (2 * radius * sin(pi / (d+1)) / (1 + x_h))^d / (d+1)^2 * L_L_evs;
end
