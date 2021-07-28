% File that will run experiments for chebyshev backward error

d = 4; % degree
polysize = 3; % size of polynomial
pkmean = 0; % mean for pseudosmith
pkwidth = 50; % std for pseudosmith

j = 1; % Chebyshev Type 1 or 2

[evs, Pmon] = polygen_pseudosmith(d, polysize, pkmean, pkwidth);

Pcheb = sym(zeros(size(Pmon)));
for r = 1:polysize
    for c = 1:polysize
        monEntry = coeffs(Pmon(r, c)); % Returns coefficients with leading coefficient last
        monEntry = fliplr(monEntry); % Returns coefficients with leading coefficient first (as mon2cheb takes it)
        chebEntryPoly = transpose(mon2cheb(monEntry, j)); % Returns coefficients with leading coefficient first IN ROW VECTOR
        chebEntryPoly = fliplr(chebEntryPoly); % Returns coefficients with leading coefficient last
        chebEntry = poly2sym(chebEntryPoly);
        Pcheb(r, c) = chebEntry;
    end
end

coeff = polygen_split_smith_mtx(d, Pcheb);

% Pcheb


