% File that will run experiments for chebyshev backward error

d = 4; % degree
polysize = 3; % size of polynomial
pkmean = 0; % mean for pseudosmith
pkwidth = 50; % std for pseudosmith

symbolic = 1;

j = 1; % Chebyshev Type 1 or 2

ep = 2;

[evs, Pmon] = polygen_pseudosmith(d, polysize, pkmean, pkwidth);

% test = polygen_split_smith_mtx(d, Pmon);
% test = double(test);
% P0 = test(:,:,1);
% P1 = test(:,:,2);
% P2 = test(:,:,3);
% P3 = test(:,:,4);
% P4 = test(:,:,5);
% polyeig(P0,P1,P2,P3,P4)
% evs

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

if symbolic == 1
    [M1sym,M0sym] = Msubfamily(d,polysize,coeff,ep,j);

    M1sym = block2notblock(M1sym);
    M0sym = block2notblock(M0sym);

    [C1sym,C0sym] = cPencil(M1sym,M0sym,j,polysize,ep,d);
    polyeig(C0sym,C1sym)
    evs
else
    coeff = double(coeff);

    [M1,M0] = Msubfamily(d,polysize,coeff,ep,j);
    M1 = block2notblock(M1);
    M0 = block2notblock(M0);

    [C1,C0] = cPencil(M1,M0,j,polysize,ep,d);

    C0;
    C1;

    polyeig(C0,C1)
    evs
    
end





