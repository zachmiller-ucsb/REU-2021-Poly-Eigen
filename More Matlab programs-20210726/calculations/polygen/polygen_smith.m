function [evs, pks] = polygen_smith(d, polysize, ...
    nodes, noderange, evgen, pkmean, pkwidth)
% Return the eigenvalues and polynomial coefficients for a matrix
% polynomial, given the size, degree, and node range.

    % obtain eigenvalues
    evs = evgen(d, nodes, noderange);

    % Build diagonal matrix for Smith form.  Note the algebraic
    % multiplicity of these eigenvalues will be 1 because their invariant
    % factors appear at the bottom right of D, consequently the eigenvalues
    % will be simple.  However, most of the matrix polynomial coefficients
    % will be singular.
    D = sym(eye(polysize));
    D(polysize, polysize) = prod(sym('x') - evs);

    % obtain random unimodular matrices
    smithleft = stdunimodular(polysize, polysize, pkmean, pkwidth);
    smithright = stdunimodular(polysize, polysize, pkmean, pkwidth);

    P = smithleft * D * smithright;
    pks = polygen_split_smith_mtx(d, P);
end
