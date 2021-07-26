function [nodes, noderange] = nodegeneval(nodegen, d)
% Evaluate the node generator with the given arguments.

    if iscell(nodegen)
        ngfunction = nodegen{1};
        ngparameters = nodegen{2};
        [nodes, noderange] = ngfunction(d, ngparameters);
        return;
    end

    [nodes, noderange] = nodegen(d);
end
