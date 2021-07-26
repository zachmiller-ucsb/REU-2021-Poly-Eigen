function answer = basiseval_monomial(d, ~, ~, evs)
% Compute the sum of the monomial basis elements evaluated at the
% eigenvalues given.  Return a row vector.

    monomial_basis = poly2sym(ones(1, d+1));
    answer = subs(monomial_basis, abs(evs));
end
