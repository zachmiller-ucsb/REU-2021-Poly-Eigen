function bounds = boundseval(boundgens, d, nodes, noderange, evs, Lx2m)
% Evaluate the bound generators with the given arguments.

    numbounds = length(boundgens);
    bounds = stdnumerize(zeros(numbounds, 2));

    for k = 1:numbounds
        [bounds(k, 1), bounds(k, 2)] = ...
            boundgens{k}(d, nodes, noderange, evs, Lx2m);
    end
end
