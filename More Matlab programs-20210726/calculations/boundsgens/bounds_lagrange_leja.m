function [u, l] = bounds_lagrange_leja(d, nodes, noderange, evs, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Use bounds derived from Leja sequences.

    x_h = max(abs(nodes));
    x_H = max(1, x_h);

    % Determine the shape of the noderange.
    [~, ~, ~, radius] = stdexplainrange(noderange);
    radius = abs(radius);

    L_H_evs = max(1 ./ (max(1, abs(evs).^d)));
    L_L_evs = min(1 ./ max(1, abs(evs)).^d);

    u = 2 * (d+1)^2 * x_H^d * L_H_evs;
    l = radius^d / (d+1) / (1+x_h)^d * L_L_evs;
end
