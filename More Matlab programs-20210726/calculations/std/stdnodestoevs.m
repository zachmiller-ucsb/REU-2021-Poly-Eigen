function [mindistances, maxdistances] = stdnodestoevs(nodes, evs)
% Return the minimum and maximum distances from the nodes to all the given
% eigenvalues.

    mindistances = zeros(1, length(evs));
    maxdistances = zeros(1, length(evs));
    for k = 1:length(evs)
        mindistances(k) = min(abs(nodes - evs(k)));
        maxdistances(k) = max(abs(nodes - evs(k)));
    end
end
