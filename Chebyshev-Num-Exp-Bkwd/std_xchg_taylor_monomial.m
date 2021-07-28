function L = std_xchg_taylor_monomial(d, ~, noderange)
% Map of bases matrix from Taylor to monomial, mapping the column
% vector of Taylor basis vectors to the column vector of monomial basis
% vectors.
%
% Note well: this is not the change of basis matrix for coordinates.

    syms x;
    [~, ~, center, ~] = stdexplainrange(noderange);
    root = -sym(center);
    L = sym(zeros(d+1));
    for k = 1:d+1
        row = sym(zeros(1, k));
        row(1, k) = sym(1);
        rootpower = root;
        for j = k-1:-1:1
            row(1:j) = row(1:j) ...
                - L(j, 1:j) * rootpower * nchoosek(sym(k-1), sym(j-1));
            rootpower = rootpower * root;
        end
        L(k, 1:k) = row;
    end
end
