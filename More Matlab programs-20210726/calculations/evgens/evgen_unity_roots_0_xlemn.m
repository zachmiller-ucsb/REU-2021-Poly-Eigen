function matrix = evgen_unity_roots_0_xlemn(n, nodes, ~)
% Return n eigenvalues consisting of the roots of unity, centered at zero,
% and scaled in radius to the minimum (nonzero) node absolute value, then
% further scaled down to zero by 10^-n.

    sortedabs = sort(abs(nodes), 'ascend');
    xl = sortedabs(1);
    if xl == 0
        % intended to fail on purpose if only one node is provided
        xl = sortedabs(2);
    end
    matrix = nodes_unity_roots_d_c_r(n, 0, xl * 10^-n);
end
