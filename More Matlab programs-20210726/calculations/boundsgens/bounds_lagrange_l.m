function [u, l] = bounds_lagrange_l(d, nodes, ~, evs, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Using Theorem 4.2 for bounds from document 2 (improvement of previous
% case for bounds_lagrange_d for lower bound and upper bound).

    x_h = max(abs(nodes));
    x_H = max(1, x_h);
    lambda_H = max(1, abs(evs));

    % Find the minimum and maximum distances between the nodes and the
    % eigenvalues.
    [~, maxnodestoevs] = stdnodestoevs(nodes, evs);

    % Find the minimum and maximum distance between the nodes. These are
    % the worst cases for the corresponding bounds.
    [minnodetonode, ~] = stdnodestonodes(nodes);

    GXd = stdgautschibound(d, nodes);

    L_H_variable = max(maxnodestoevs.^d ./ (1 + abs(evs).^d));
    L_H_choice = (d+1) / minnodetonode^d * L_H_variable;
    L_L_choice = (1 / (d+1)) * min(1 ./ lambda_H.^d);

    u = (d+1) * x_H^d * L_H_choice;
    l = L_L_choice / ((d+1) * GXd);
end
