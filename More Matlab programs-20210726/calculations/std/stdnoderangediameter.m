function diameter = stdnoderangediameter(noderange)
% Return the diameter of the ball enclosing the given noderange.
% This method assumes the node generators are kind enough to facilitate.

    absvalues = abs(noderange);
    [~, sortidx] = sort(absvalues, 'descend');
    diameter = abs(noderange(sortidx(1)) ...
        - noderange(sortidx(length(noderange))));
end
