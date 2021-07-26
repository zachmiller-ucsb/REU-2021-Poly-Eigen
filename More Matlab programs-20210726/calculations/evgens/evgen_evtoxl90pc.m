function matrix = evgen_evtoxl90pc(n, nodes, noderange)
% Return n equidistant eigenvalues in the line between the center of the
% noderange and the node with largest absolute value.

    [~, ~, center, ~] = stdexplainrange(noderange);
    absnodes = abs(nodes);
    [~, sortedabsindices] = sort(absnodes, 'ascend');
    xl = nodes(sortedabsindices(1));
    matrix = xl + sqrt((1:n) / (10*n)) * (center - xl);
end
