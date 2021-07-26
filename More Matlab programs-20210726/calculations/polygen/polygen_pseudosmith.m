function [evs, pks] = polygen_pseudosmith(d, polysize, ...
    nodes, noderange, evgen, pkmean, pkwidth)
% Return the eigenvalues and polynomial coefficients for a matrix
% polynomial, given the size, degree, and node range.

    % This polynomial generator builds an upper triangular matrix UT with
    % random polynomials above the diagonal, and characteristic polynomial
    % factors of the form (x - ev) in the diagonal.  Then, the UT matrix is
    % mixed by two random scalar unimodular matrices.  Both to prevent
    % infinite eigenvalues as well as to produce realistic matrix
    % polynomials, the matrix polynomial will have d * polysize factors in
    % the diagonal.  When this is done, trying to reduce the scalar
    % matrices corresponding to each power of the independent variable will
    % (virtually always) result in the identity matrix.

    % obtain eigenvalues
    evs = evgen(d * polysize, nodes, noderange);

    % Start from something that looks like the diagonal matrix in Smith's
    % normal form.  However, this process builds an upper triangular
    % matrix instead to prevent the matrix polynomial coefficients from
    % being singular when they correspond to degrees larger than zero.
    UT = sym(zeros(polysize));
    for k = 1:polysize
        evrowindices = (k-1) * d + 1 : k * d;
        UT(k, k) = prod(sym('x') - evs(evrowindices));
    end

    % Deviate from Smith further here to avoid singular coefficients.
    % Fill the upper triangular portion of UT with random polynomials.
    for r = 1:polysize-1
        for c = r+1:polysize
            poly = poly2sym(stdrandn(1, d+1, pkwidth) + pkmean);
            while length(coeffs(poly)) ~= d+1
                poly = poly2sym(stdrandn(1, d+1, pkwidth) + pkmean);
            end
            UT(r, c) = poly;
        end
    end

    % Obtain random unimodular matrices.  Use normals between 0 and 1, the
    % matrix UT already has entry coefficients with the requested normal
    % distribution.  Note this does not disturb the eigenvalues because
    % det(SL * UT * SR) = det(SL) * det(UT) * det(SR), and the outer
    % determinants are constants because SL and SR are scalar matrices.
    smithleft = stdunimodular(polysize, polysize, 0, 1);
    smithright = stdunimodular(polysize, polysize, 0, 1);

    P = smithleft * UT * smithright;
    pks = polygen_split_smith_mtx(d, P);
end
