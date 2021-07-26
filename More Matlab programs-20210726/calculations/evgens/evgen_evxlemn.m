function matrix = evgen_evxlemn(n, nodes, ~)
% Return n eigenvalues consisting of the smallest non-zero node absolute
% value, then further scaled down to zero by 10^-1 through 10^-n.

    sortedabs = sort(abs(nodes), 'ascend');
    xl = sortedabs(1);
    if xl == 0
        % intended to fail on purpose if only one node is provided
        xl = sortedabs(2);
    end
    matrix = xl * sym(10).^(-n:-1);
end
