function answer = basiseval_chebyshev(d, ~, ~, evs)
% Compute the sum of the Chebyshev basis elements evaluated at the
% eigenvalues given.  Return a row vector.

    cps = std_chebyshev_polynomials(d);
    answer = sym(zeros(size(evs)));
    for k = 1:d+1
        cpsym = poly2sym(cps{k});
        each = subs(cpsym, evs);
        answer = answer + abs(each);
    end
end
