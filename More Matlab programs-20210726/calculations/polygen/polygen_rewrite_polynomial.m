function polyout = polygen_rewrite_polynomial(L, polyin)
% This function takes a set of matrix coefficients for a polynomial written
% in basis A, and outputs a new set of matrix coefficients for the
% polynomial written in basis B, by applying the basis change matrix L over
% the implied variable parts of each monomial in the matrix polynomial.
%
% It is implied that the polynomial given has degree d, that L has
% dimension (d+1)x(d+1), and that hence there are d+1 matrix coefficients
% to rewrite using L.  Further, it is assumed polyin is a 3d array of size
% mxmx(d+1), with the polynomial coefficients indexed by (:,:,k) with
% 1 <= k <= d+1, and m being the size the (square) given matrix polynomial.
%
% Importantly, L is expected to be the change of basis matrix that maps
% bases in the opposite direction as needed for these coefficients.  That
% is, L is the change of basis matrix that takes the basis vectors for the
% basis B to the basis vectors of basis A.  Also, the matrix L must be
% constructed and provided here in symbolic form.

    polyout = sym(zeros(size(polyin)));
    polylength = size(polyin, 3);
    rewritten = cell(polylength);

    for k = 1:polylength
        polyk = polyin(:,:,k);
        for j = 1:polylength
            if k == 1
                rewritten{j} = L(k, j) * polyk;
            else
                rewritten{j} = rewritten{j} + L(k, j) * polyk;
            end
        end
    end
    for j = 1:polylength
        polyout(:,:,j) = rewritten{j};
    end
end
