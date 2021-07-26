function [u, l] = bounds_lagrange_d(d, nodes, ~, evs, ~)
% Calculate a Lagrange specific bound for the ratio of condition numbers.
% Using Theorem 4.2 for bounds from document 2.
% Using Lemma 21 for L_L and L_H expressions, choosing the branch that
% raises to the power of d.

    x_h = max(abs(nodes));
    x_H = max(1, x_h);

    % Find the minimum and maximum distances between the nodes and the
    % eigenvalues.
    [minnodestoevs, maxnodestoevs] = stdnodestoevs(nodes, evs);

    % Find the minimum and maximum distance between the nodes. These are
    % the worst cases for the corresponding bounds.
    [minnodetonode, maxnodetonode] = stdnodestonodes(nodes);

    GXd = stdgautschibound(d, nodes);

    L_H_variable = max(maxnodestoevs.^d ./ (1 + abs(evs).^d));
    L_H_choice = (d+1) / minnodetonode^d * L_H_variable;
    L_L_variable = min(minnodestoevs.^d ./ (max(1, abs(evs))).^d);
    L_L_choice = 1 / maxnodetonode^d * L_L_variable;

    u = (d+1) * x_H^d * L_H_choice;
    l = L_L_choice / ((d+1)^2 * GXd);
end
