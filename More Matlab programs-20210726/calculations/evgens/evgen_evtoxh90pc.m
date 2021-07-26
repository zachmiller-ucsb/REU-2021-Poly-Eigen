function matrix = evgen_evtoxh90pc(n, nodes, noderange)
% Return n equidistant eigenvalues in the line between the center of the
% noderange and the node with largest absolute value.

    [~, ~, center, ~] = stdexplainrange(noderange);
    absnodes = abs(nodes);
    [~, sortedabsindices] = sort(absnodes, 'descend');
    xh = nodes(sortedabsindices(1));
    matrix = center + sqrt((9*n+1:10*n) / (10*n+1)) * (xh - center);
end
