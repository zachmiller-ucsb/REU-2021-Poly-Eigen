function L = std_xchg_lagrange_monomial(~, nodes, ~)
% Map of bases matrix from Lagrange to monomial, mapping the column
% vector of Lagrange basis vectors to the column vector of monomial basis
% vectors.
%
% Note well: this is not the change of basis matrix for coordinates.

    L = stdvander(sym(nodes)).';
end
