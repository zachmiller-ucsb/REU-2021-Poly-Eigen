function cp = std_chebyshev_polynomials(d)
% Return the Chebyshev polynomials of the first kind, with degrees from 0
% to d, as a cell array of coefficients in the monomial basis.

    % use a cell array because the polynomial coefficient matrices have
    % different sizes
    cp = cell(d+1, 1);
    cp{1} = 1;      % T_0 = 1
    cp{2} = [1 0];  % T_1 = x
    for k = 3:d+1   % T_n = 2xT_n-1 - T_n-2
        cp{k} = [2*cp{k-1} 0] - [0 0 cp{k-2}];
    end
end
