function [u, l] = bounds_newton(d, nodes, ~, evs, ~)
% Calculate a Newton specific bound for the ratio of condition numbers.

    xh = max(abs(nodes));

    % Find the minimum and maximum distances between the nodes and the
    % eigenvalues.
    [minnodestoevs, maxnodestoevs] = stdnodestoevs(nodes, evs);

    evuppernumerators = max(1, maxnodestoevs.^d);
    evupperdenominators = 1 + abs(evs).^d;
    evupperfractions = evuppernumerators ./ evupperdenominators;

    evlowernumerators = max(1, minnodestoevs.^d);
    evlowerdenominators = max(1, abs(evs).^d);
    evlowerfractions = evlowernumerators ./ evlowerdenominators;

    u = (d+1)^2 * (1+xh)^d * max(evupperfractions);
    l = min(evlowerfractions) / (d+1)^2 / (1+xh)^d;
end
