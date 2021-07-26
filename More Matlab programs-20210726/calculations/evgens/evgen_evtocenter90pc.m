function matrix = evgen_evtocenter90pc(n, nodes, noderange)
% Return n equidistant eigenvalues in the line between the center of the
% noderange and the node with largest absolute value.

    [~, ~, center, ~] = stdexplainrange(noderange);
    absnodes = abs(nodes);
    [~, sortedabsindices] = sort(absnodes, 'descend');
    xh = nodes(sortedabsindices(1));
    matrix = center + ((1:n) / (10*n)).^(3/2) * (xh - center);
end
