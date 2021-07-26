function matrix = evgen_evtoeq(n, nodes, noderange)
% Return n equidistant eigenvalues in the line between the center of the
% noderange and the node with largest absolute value.

    [~, ~, center, ~] = stdexplainrange(noderange);
    absnodes = abs(nodes);
    [~, sortedabsindices] = sort(absnodes, 'descend');
    xh = nodes(sortedabsindices(1));
    matrix = center + (1:n) / (n+1) * (xh - center);
end
