function matrix = evgen_evdipole99pc(n, nodes, ~)
% Return n eigenvalues in the line between the node with largest absolute
% value and the node furthest away from that first node, bunched up towards
% both nodes.

    absnodes = abs(nodes);
    [~, sortedabsindices] = sort(absnodes, 'descend');
    xh_index = sortedabsindices(1);
    xh = nodes(xh_index);
    nodes_noxh = nodes([1:xh_index-1 xh_index+1:length(nodes)]);
    distances = abs(nodes_noxh - xh);
    [~, sorteddistancesindices] = sort(distances, 'descend');
    xl = nodes_noxh(sorteddistancesindices(1));
    nl = floor(n/2);
    nh = ceil(n/2);
    mh = xl + sqrt((199*nh+1:200*nh) / (200*nh+1)) * (xh - xl);
    ml = fliplr(mh);
    ml = xh + xl - ml(1:nl);
    matrix = [ml mh];
end
