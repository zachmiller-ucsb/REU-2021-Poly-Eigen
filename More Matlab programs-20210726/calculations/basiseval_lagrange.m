function answer = basiseval_lagrange(d, nodes, ~, evs)
% Compute the sum of the Lagrange basis elements evaluated at the
% eigenvalues given.  Return a row vector.

    syms x;
    symnodes = sym(nodes);
    answer = sym(zeros(size(evs)));
    for k = 1:d+1
        xk = symnodes(k);
        nodes_noxk = nodes([1:k-1 k+1:d+1]);
        lknumerator = prod(x - nodes_noxk, 2);
        lkdenominator = prod(xk - nodes_noxk, 2);
        answer = answer + abs(subs(lknumerator, evs) / lkdenominator);
    end
end