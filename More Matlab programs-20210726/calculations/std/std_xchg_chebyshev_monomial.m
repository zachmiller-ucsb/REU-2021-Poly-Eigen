function L = std_xchg_chebyshev_monomial(d, j, ~, ~)
% Map of bases matrix from Chebyshev to monomial, mapping the column
% vector of Chebyshev basis vectors to the column vector of monomial basis
% vectors.
%
% Note well: this is not the change of basis matrix for coordinates.

    cps = std_chebyshev_polynomials(d,j);
    L = zeros(d+1);
    for k = 1:d+1
        L(k:d+1, d+2-k) = cps{d+2-k};
    end
    L = flipud(L);
    L = inv(sym(L.'));
end
