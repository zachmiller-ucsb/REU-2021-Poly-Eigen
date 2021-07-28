function L = std_xchg_newton_monomial(d, nodes, ~)
% Map of basis matrix from Newton to monomial, mapping the column
% vector of Newton basis vectors to the column vector of monomial basis
% vectors.
%
% Calculate L directly via Gander's method.  This is much faster than
% either the direct expression of H functions, or factoring V = UL then
% inverting some relevant matrix.
%
% Note well: this is not the change of basis matrix for coordinates.

    syms x;
    symnodes = sym(nodes);
    V = stdvander(symnodes);

    invU = sym(zeros(d+1));
    for c = 1:d+1
        for r = 1:c
            crnodes = symnodes([1:r-1 r+1:c]);
            npoly_xr = prod(symnodes(r) - crnodes, 2);
            invU(r, c) = 1 / npoly_xr;
        end
    end
    L = V.' * invU;  % i.e., L = V^T * U^-1
end
