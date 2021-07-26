function [mindistance, maxdistance] = stdnodestonodes(nodes)
% Return the minimum and maximum distance between the given nodes.

    maxdistance = abs(nodes(2) - nodes(1));
    mindistance = maxdistance;
    for k = 1:length(nodes)
       nodes_noxk = nodes([1:k-1 k+1:length(nodes)]);
       distances = abs(nodes_noxk - nodes(k));
       maxdistance = max(maxdistance, max(distances));
       mindistance = min(mindistance, min(distances));
    end
end
