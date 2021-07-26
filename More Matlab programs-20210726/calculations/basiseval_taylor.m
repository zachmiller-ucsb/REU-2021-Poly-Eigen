function answer = basiseval_taylor(d, ~, noderange, evs)
% Compute the sum of the Taylor basis elements evaluated at the
% eigenvalues given (i.e. the monomial basis shifted to the geometric
% center of the nodes).  Return a row vector.

    syms x;
    [~, ~, center, ~] = stdexplainrange(noderange);
    factor = x - center;
    basis_polynomial = sym(1);
    answer = stdnumerize(ones(1, length(evs)));
    for k = 1:d
        basis_polynomial = basis_polynomial * factor;
        answer = answer + abs(subs(basis_polynomial, evs));
    end
end
