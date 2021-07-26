function answer = basiseval_newton(d, nodes, ~, evs)
% Compute the sum of the Newton basis elements evaluated at the
% eigenvalues given.  Return a row vector.

    syms x;
    symnodes = sym(nodes);
    answer = sym(zeros(size(evs)));
    nk = sym(1);
    for k = 1:d+1
        answer = answer + abs(subs(nk, evs));
        nk = nk * (x - symnodes(k));
    end
end
