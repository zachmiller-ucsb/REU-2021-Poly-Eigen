function pks = polygen_split_smith_mtx(d, P)
% Given the output of the Smith process to generate matrix polynomials,
% return the matrix coefficients for the matrix polynomial.  The last
% coefficient is the matrix polynomial's head coefficient.

    n = size(P, 1);
    pks = sym(zeros(n, n, d+1));
    for r = 1:n
        for c = 1:n
            pks(r, c, :) = coeffs(P(r, c));
        end
    end
end
